#!/usr/bin/env python3
"""Repair exons after miniprot.

1. Leader: a short GT-AG exon that already contains methionine. Do not prepend
   intron sequence, and do not walk an intron looking for an ATG.
2. I exon: tblastn of the I-domain peptides, snapped to AG…GT and spliced in.
3. An I hit with no gene is kept as an I-only relic.

A BLAST alignment is not an exon. Exon ends are splice sites.
"""
from __future__ import annotations

import subprocess
from collections import defaultdict
from pathlib import Path

import domain_annotation as ann
import gene_models as p

HYDRO = set("AILMFWV")
LEAD_INTRON = (40, 2500)  # leader-to-V introns are usually a few hundred bp
LEAD_NT = (39, 90)  # observed leader exons are 43–58 nt (phase 1 into V)
I_FLANK = 4000
MIN_I_AA = 70
MAX_I_AA = 130
I_NT = (210, 400)
SNAP = 25
I_SNAP = 60
V_AG_EXTEND = 300


def _gene(row: dict) -> dict:
    exons = list(row.get("exons") or [])
    return {
        "contig": row["contig"],
        "strand": row.get("strand") or "+",
        "cds": [{"start": int(e["start"]), "end": int(e["end"]), "phase": str(e.get("phase") or "0")} for e in exons],
    }


def rebuild(row: dict, genome: dict[str, str], note: str | None = None) -> dict:
    gene = _gene(row)
    nt, exons, notes = p.spliced_cds(genome, gene)
    aa = p.translate(nt) if nt else ""
    stops_mid = "*" in aa[:-1] if aa else True
    row = dict(row)
    row["exons"] = exons
    row["aa"] = aa
    row["nt_len"] = len(nt)
    row["splice_ok"] = bool(nt) and not any(n.startswith("noncanonical") for n in notes)
    row["orf_ok"] = bool(aa) and (not stops_mid) and len(aa) >= 40
    row["notes"] = list(row.get("notes") or [])
    if isinstance(row["notes"], str):
        row["notes"] = [x for x in row["notes"].split(";") if x and x != "."]
    if note:
        row["notes"].append(note)
    row["notes"].extend(notes)
    if exons:
        row["start"] = min(int(e["start"]) for e in exons)
        row["end"] = max(int(e["end"]) for e in exons)
        row["span"] = int(row["end"]) - int(row["start"]) + 1
        row["n_exons"] = len(exons)
    return row


def _overlap(a1: int, a2: int, b1: int, b2: int) -> int:
    return max(0, min(a2, b2) - max(a1, b1) + 1)


def tblastn(query: Path, genome_fa: Path, dest: Path, evalue: str = "1e-5") -> list[dict]:
    dest.parent.mkdir(parents=True, exist_ok=True)
    db = genome_fa.with_suffix(".blastdb")
    if not list(db.parent.glob(db.name + ".n*")):
        subprocess.check_call(
            ["makeblastdb", "-in", str(genome_fa), "-dbtype", "nucl", "-out", str(db)],
            stdout=subprocess.DEVNULL,
        )
    with dest.open("w") as fh:
        subprocess.check_call(
            [
                "tblastn",
                "-query",
                str(query),
                "-db",
                str(db),
                "-evalue",
                evalue,
                "-seg",
                "no",
                "-num_threads",
                p.THREADS,
                "-outfmt",
                "6 qseqid sseqid pident length qlen sstart send evalue bitscore sstrand",
            ],
            stdout=fh,
        )
    hits = []
    if not dest.exists():
        return hits
    for line in dest.read_text().splitlines():
        col = line.split("\t")
        if len(col) < 10:
            continue
        s, e = int(col[5]), int(col[6])
        hits.append(
            {
                "q": col[0],
                "contig": col[1],
                "pident": float(col[2]),
                "length": int(col[3]),
                "qlen": int(col[4]),
                "start": min(s, e),
                "end": max(s, e),
                "evalue": float(col[7]),
                "bits": float(col[8]),
                "strand": "+" if col[9] in {"plus", "+"} else "-",
            }
        )
    return hits


def merge_hits(hits: list[dict], slack: int = 30) -> list[dict]:
    by = defaultdict(list)
    for h in hits:
        by[(h["contig"], h["strand"])].append(h)
    out = []
    for (ctg, strand), grp in by.items():
        grp = sorted(grp, key=lambda x: (x["start"], -x["bits"]))
        cur = None
        for h in grp:
            if cur is None or h["start"] > cur["end"] + slack:
                if cur:
                    out.append(cur)
                cur = dict(h)
            else:
                cur["end"] = max(cur["end"], h["end"])
                cur["bits"] = max(cur["bits"], h["bits"])
                cur["pident"] = max(cur["pident"], h["pident"])
        if cur:
            out.append(cur)
    return out


def _has_i_exon(row: dict) -> bool:
    """True if a middle exon is already I-sized. HMM n_i is not enough."""
    exons = list(row.get("exons") or [])
    if len(exons) < 2:
        return False
    for e in exons[1:]:
        aa_len = (int(e["end"]) - int(e["start"]) + 1) // 3
        if MIN_I_AA <= aa_len <= MAX_I_AA:
            return True
    return False


def _is_signal(aa: str) -> bool:
    pep = (aa or "").replace("*", "")
    if not pep.startswith("M") or not (12 <= len(pep) <= 28):
        return False
    if sum(c in HYDRO for c in pep[1:16]) < 6:
        return False
    # First 22 aa of a V domain is not a signal peptide.
    if "GSGSR" in pep or "EGTSL" in pep or "EGYIL" in pep:
        return False
    return True


def _pident(a: str, b: str) -> float:
    a = (a or "")[:20]
    b = (b or "")[:20]
    best = 0.0
    for i in range(len(a)):
        for j in range(len(b)):
            n = min(len(a) - i, len(b) - j)
            if n < 8:
                break
            ident = sum(a[i + k] == b[j + k] for k in range(n)) / n
            if ident > best:
                best = ident
    return best


def _leader_score(aa: str, known: list[str]) -> float:
    pep = (aa or "").replace("*", "")
    if not known:
        return 0.0
    return max(_pident(pep, k) for k in known)


def _v_acceptor_ok(seq: str, first: dict, strand: str) -> bool:
    fs, fe = int(first["start"]), int(first["end"])
    if strand == "+":
        return seq[fs - 3 : fs - 1].upper() == "AG"
    return p.revcomp(seq[fe : fe + 2]).upper() == "AG"


def snap_v_acceptor(row: dict, genome: dict[str, str]) -> dict:
    """Move the 5' end of exon 1 to the nearest AG. miniprot often starts mid-V."""
    exons = list(row.get("exons") or [])
    contig = row.get("contig")
    if not exons or contig not in genome:
        return row
    seq = genome[contig]
    strand = row.get("strand") or "+"
    first = dict(exons[0])
    fs, fe = int(first["start"]), int(first["end"])
    if _v_acceptor_ok(seq, first, strand):
        return row
    aa0 = (row.get("aa") or "").replace("*", "")
    best = None
    if strand == "+":
        for new_fs in range(fs, max(1, fs - V_AG_EXTEND), -1):
            if seq[new_fs - 3 : new_fs - 1].upper() != "AG":
                continue
            trial = dict(row)
            ex0 = dict(first)
            ex0["start"] = new_fs
            trial["exons"] = [ex0] + exons[1:]
            trial = rebuild(trial, genome)
            new_aa = (trial.get("aa") or "").replace("*", "")
            if trial.get("splice_ok") and trial.get("orf_ok") and len(new_aa) >= max(40, len(aa0) - 5):
                best = trial
                break
    else:
        n = len(seq)
        for new_fe in range(fe, min(n - 2, fe + V_AG_EXTEND) + 1):
            if p.revcomp(seq[new_fe : new_fe + 2]).upper() != "AG":
                continue
            trial = dict(row)
            ex0 = dict(first)
            ex0["end"] = new_fe
            trial["exons"] = [ex0] + exons[1:]
            trial = rebuild(trial, genome)
            new_aa = (trial.get("aa") or "").replace("*", "")
            if trial.get("splice_ok") and trial.get("orf_ok") and len(new_aa) >= max(40, len(aa0) - 5):
                best = trial
                break
    if best is None:
        return row
    best["notes"] = list(best.get("notes") or []) + ["V_acceptor_snapped"]
    return best


def _leader_hsp_ok(h: dict, first: dict, strand: str) -> bool:
    fs, fe = int(first["start"]), int(first["end"])
    lo, hi = LEAD_INTRON
    if strand == "+":
        if h["end"] >= fs:
            return False
        gap = fs - h["end"]
        return lo <= gap <= hi
    if h["start"] <= fe:
        return False
    gap = h["start"] - fe
    return lo <= gap <= hi


def _snap_leader_plus(seq: str, first: dict, hsp: dict) -> list[dict]:
    fs = int(first["start"])
    n = len(seq)
    if seq[fs - 3 : fs - 1].upper() != "AG":
        return []
    win_s = max(1, hsp["start"] - SNAP)
    win_e = min(hsp["end"] + SNAP, fs - LEAD_INTRON[0])
    out = []
    lo, hi = LEAD_NT
    for pos in range(win_s, max(win_s, win_e - lo + 1) + 1):
        if seq[pos - 1 : pos + 2].upper() != "ATG":
            continue
        for end in range(pos + lo - 1, min(pos + hi - 1, win_e) + 1):
            if seq[end : end + 2].upper() != "GT":
                continue
            intron = fs - end - 1
            if not (LEAD_INTRON[0] <= intron <= LEAD_INTRON[1]):
                continue
            aa = p.translate(seq[pos - 1 : end])
            if not _is_signal(aa):
                continue
            out.append({"start": pos, "end": end, "phase": "0", "hydro": sum(c in HYDRO for c in aa[1:16]), "aa": aa})
    return out


def _snap_leader_minus(seq: str, first: dict, hsp: dict) -> list[dict]:
    fe = int(first["end"])
    n = len(seq)
    if p.revcomp(seq[fe : fe + 2]).upper() != "AG":
        return []
    win_s = max(fe + LEAD_INTRON[0], hsp["start"] - SNAP)
    win_e = min(n, hsp["end"] + SNAP)
    out = []
    lo, hi = LEAD_NT
    for le in range(max(win_s + lo - 1, fe + LEAD_INTRON[0] + lo), win_e + 1):
        if p.revcomp(seq[le - 3 : le]).upper() != "ATG":
            continue
        for ls in range(le - hi + 1, le - lo + 2):
            if ls <= fe + LEAD_INTRON[0]:
                continue
            if p.revcomp(seq[ls - 3 : ls - 1] if ls >= 3 else "").upper() != "GT":
                continue
            intron = ls - fe - 1
            if not (LEAD_INTRON[0] <= intron <= LEAD_INTRON[1]):
                continue
            aa = p.translate(p.revcomp(seq[ls - 1 : le]))
            if not _is_signal(aa):
                continue
            out.append({"start": ls, "end": le, "phase": "0", "hydro": sum(c in HYDRO for c in aa[1:16]), "aa": aa})
    return out


def snap_leader_exon(genome: dict[str, str], contig: str, strand: str, first: dict, hsp: dict) -> dict | None:
    seq = genome.get(contig, "")
    if not seq:
        return None
    cands = _snap_leader_plus(seq, first, hsp) if strand == "+" else _snap_leader_minus(seq, first, hsp)
    if not cands:
        return None
    cands.sort(key=lambda x: (-x["hydro"], -(x["end"] - x["start"])))
    best = cands[0]
    return {"start": best["start"], "end": best["end"], "phase": "0"}


def _snap_i_plus(seq: str, hsp: dict) -> list[dict]:
    n = len(seq)
    hs, he = hsp["start"], hsp["end"]
    out = []
    lo, hi = I_NT
    for start in range(max(3, hs - I_SNAP), min(n, hs + I_SNAP) + 1):
        if seq[start - 3 : start - 1].upper() != "AG":
            continue
        for end in range(max(start + lo - 1, he - I_SNAP), min(n - 2, he + I_SNAP, start + hi - 1) + 1):
            if seq[end : end + 2].upper() != "GT":
                continue
            out.append({"start": start, "end": end, "phase": "0"})
    return out


def _snap_i_minus(seq: str, hsp: dict) -> list[dict]:
    n = len(seq)
    hs, he = hsp["start"], hsp["end"]
    out = []
    lo, hi = I_NT
    for end in range(max(1, he - I_SNAP), min(n - 2, he + I_SNAP) + 1):
        if p.revcomp(seq[end : end + 2]).upper() != "AG":
            continue
        for start in range(max(3, hs - I_SNAP, end - hi + 1), min(hs + I_SNAP, end - lo + 1) + 1):
            if p.revcomp(seq[start - 3 : start - 1]).upper() != "GT":
                continue
            out.append({"start": start, "end": end, "phase": "0"})
    return out


def snap_i_exon(genome: dict[str, str], contig: str, strand: str, hsp: dict) -> list[dict]:
    seq = genome.get(contig, "")
    if not seq:
        return []
    return _snap_i_plus(seq, hsp) if strand == "+" else _snap_i_minus(seq, hsp)


def recover_leader(row: dict, genome: dict[str, str], lead_hits: list[dict], known: list[str]) -> dict:
    aa = (row.get("aa") or "").replace("*", "")
    exons = list(row.get("exons") or [])
    if not exons or aa.startswith("M") or row["contig"] not in genome:
        return row
    if int(row.get("n_v") or 0) < 1:
        return row
    strand = row.get("strand") or "+"
    first = exons[0]
    contig = row["contig"]
    seq = genome[contig]
    fs, fe = int(first["start"]), int(first["end"])
    windows = [h for h in lead_hits if h["contig"] == contig and h["strand"] == strand and _leader_hsp_ok(h, first, strand)]
    if strand == "+":
        windows.append({"start": max(1, fs - LEAD_INTRON[1]), "end": max(1, fs - LEAD_INTRON[0]), "bits": 0})
    else:
        windows.append({"start": fe + LEAD_INTRON[0], "end": min(len(seq), fe + LEAD_INTRON[1]), "bits": 0})
    snap_fn = _snap_leader_plus if strand == "+" else _snap_leader_minus
    seen: set[tuple[int, int]] = set()
    best = None
    best_sc = -1.0
    for h in windows:
        for new_ex in snap_fn(seq, first, h):
            key = (int(new_ex["start"]), int(new_ex["end"]))
            if key in seen or _overlap(key[0], key[1], fs, fe) >= 10:
                continue
            seen.add(key)
            sc = _leader_score(new_ex.get("aa") or "", known)
            if sc < 0.30:
                continue
            trial = dict(row)
            trial["exons"] = [new_ex] + exons
            trial = rebuild(trial, genome)
            new_aa = (trial.get("aa") or "").replace("*", "")
            if not trial.get("splice_ok"):
                continue
            if not new_aa.startswith("M"):
                ann.apply_leader_trim(trial)
                new_aa = (trial.get("aa") or "").replace("*", "")
            if not (new_aa.startswith("M") and trial.get("orf_ok") and len(new_aa) >= len(aa)):
                continue
            if sc > best_sc:
                trial["notes"] = list(trial.get("notes") or []) + ["leader_exon_recovered"]
                best = trial
                best_sc = sc
    if best is not None and best_sc >= 0.30:
        return best
    return row


def recover_i_exon(row: dict, genome: dict[str, str], i_hits: list[dict]) -> dict:
    exons = list(row.get("exons") or [])
    if not exons or _has_i_exon(row) or row["contig"] not in genome:
        return row
    strand = row.get("strand") or "+"
    contig = row["contig"]
    span_s, span_e = int(row["start"]), int(row["end"])
    cands = []
    for h in i_hits:
        if h["contig"] != contig or h["strand"] != strand:
            continue
        aa_len = (h["end"] - h["start"] + 1) // 3
        if aa_len < MIN_I_AA - 10 or aa_len > MAX_I_AA + 20:
            continue
        if _overlap(h["start"], h["end"], span_s - I_FLANK, span_e + I_FLANK) == 0:
            continue
        if any(_overlap(h["start"], h["end"], int(e["start"]), int(e["end"])) >= 40 for e in exons):
            continue
        cands.append(h)
    cands.sort(key=lambda x: -x["bits"])
    for h in cands:
        for new_ex in snap_i_exon(genome, contig, strand, h):
            if any(_overlap(new_ex["start"], new_ex["end"], int(e["start"]), int(e["end"])) >= 40 for e in exons):
                continue
            cds = sorted(exons + [new_ex], key=lambda e: int(e["start"]), reverse=(strand == "-"))
            trial = dict(row)
            trial["exons"] = cds
            trial = rebuild(trial, genome)
            if not trial.get("splice_ok") or not trial.get("orf_ok"):
                continue
            trial["notes"] = list(trial.get("notes") or []) + ["I_exon_recovered"]
            return trial
    return row


def orphan_i_relics(unique: list[dict], i_hits: list[dict], genome: dict[str, str], gid: str) -> list[dict]:
    relics = []
    n = 0
    for h in i_hits:
        aa_len = (h["end"] - h["start"] + 1) // 3
        if aa_len < MIN_I_AA - 15 or aa_len > MAX_I_AA + 25:
            continue
        if h["contig"] not in genome:
            continue
        covered = False
        for r in unique:
            if r["contig"] != h["contig"]:
                continue
            if _overlap(h["start"], h["end"], int(r["start"]), int(r["end"])) >= 40:
                covered = True
                break
        if covered:
            continue
        n += 1
        row = {
            "gene_id": f"{gid}_i_{n:04d}",
            "source": "i_tblastn_relic",
            "contig": h["contig"],
            "start": h["start"],
            "end": h["end"],
            "strand": h["strand"],
            "exons": [{"start": h["start"], "end": h["end"], "phase": "0"}],
            "query": h.get("q", "."),
            "identity": f"{h.get('pident', 0):.3f}",
            "notes": ["orphan_I"],
            "look_alike": ".",
            "known": ".",
            "gff_hint": ".",
        }
        row = rebuild(row, genome, "orphan_I")
        relics.append(row)
    print(f"Orphan I relics added: {len(relics)}", flush=True)
    return relics


def write_leader_fa(proteins: list[tuple[str, str]], dest: Path) -> list[str]:
    recs = []
    seen = set()
    known: list[str] = []
    for name, aa in proteins:
        pep = (aa or "").replace("*", "")
        picked = ""
        for L in (14, 16, 18, 20, 22):
            if _is_signal(pep[:L]):
                picked = pep[:L]
                break
        if not picked or picked in seen:
            continue
        seen.add(picked)
        recs.append((str(name) + "_lead", picked))
        known.append(picked)
    p.write_fasta(dest, recs)
    return known


def repair_loci(
    unique: list[dict],
    genome: dict[str, str],
    genome_fa: Path,
    i_fa: Path,
    known_arlee: Path,
    out: Path,
    gid: str,
) -> list[dict]:
    lead_fa = out / "seeds" / "leader40.faa"
    lead_src: list[tuple[str, str]] = []
    if known_arlee.exists():
        for hdr, aa in p.read_fasta(known_arlee).items():
            lead_src.append((p.fasta_id(hdr), aa))
    for row in unique:
        if row.get("status") != "complete":
            continue
        exons = list(row.get("exons") or [])
        if not exons or row["contig"] not in genome:
            continue
        e0 = exons[0]
        nt = p.slice_seq(genome, row["contig"], int(e0["start"]), int(e0["end"]))
        if (row.get("strand") or "+") == "-":
            nt = p.revcomp(nt)
        pep = p.translate(nt).replace("*", "")
        if _is_signal(pep):
            lead_src.append((row["gene_id"] + "_ex1", pep))
    known = write_leader_fa(lead_src, lead_fa)
    print(f"Leader peptides for homology: {len(known)}", flush=True)
    lead_hits = []
    if lead_fa.exists() and lead_fa.stat().st_size > 10:
        print("+ tblastn leader peptides", flush=True)
        lead_hits = tblastn(lead_fa, genome_fa, out / "miniprot" / "leader_vs_genome.tblastn", evalue="1e-3")
        lead_hits = [h for h in lead_hits if 8 <= h["length"] <= 28]
        print(f"Leader HSPs: {len(lead_hits)}", flush=True)
    print("+ tblastn published I domains", flush=True)
    i_raw = tblastn(i_fa, genome_fa, out / "miniprot" / "I_vs_genome.tblastn", evalue="1e-5")
    i_hits = merge_hits([h for h in i_raw if h["length"] >= 60])
    print(f"I intervals: {len(i_hits)} (raw HSPs {len(i_raw)})", flush=True)

    n_v = n_lead = n_i = 0
    fixed = []
    for row in unique:
        a = snap_v_acceptor(row, genome)
        if a is not row and "V_acceptor_snapped" in (a.get("notes") or []):
            n_v += 1
        b = recover_leader(a, genome, lead_hits, known)
        if b is not a and "leader_exon_recovered" in (b.get("notes") or []):
            n_lead += 1
        c = recover_i_exon(b, genome, i_hits)
        if c is not b and "I_exon_recovered" in (c.get("notes") or []):
            n_i += 1
        fixed.append(c)
    print(f"Snapped V acceptors={n_v} Recovered leaders={n_lead} I-exons={n_i}", flush=True)
    relics = orphan_i_relics(fixed, i_hits, genome, gid)
    return fixed + relics
