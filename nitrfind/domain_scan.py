#!/usr/bin/env python3
"""Domain scan of a six-frame translation.

Index open reading frames from esl-translate, turn HMMER hits into genomic
coordinates, and cut a window around each locus that does not already have a
miniprot model.
"""
from __future__ import annotations

import re
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import gene_models as p

THREADS = p.THREADS

ORF_RE = re.compile(
    r"^>(orf\d+)\s+source=(\S+)\s+coords=(\d+)\.\.(\d+)\s+length=(\d+)\s+frame=(\d+)"
)
FLANK = 4000
MAX_VI_GAP = 4000
SCORE_MIN = p.SCORE_MIN
SPAN_CAP = p.SPAN_CAP
MAX_ORF_AA = 2000


def filter_orfs(src: Path, dest: Path, max_aa: int = MAX_ORF_AA) -> None:
    """hmmsearch aborts if any target is >100K aa. NITR domains are ~100 aa."""
    kept = dropped = 0
    keep = False
    print(f"Filtering ORFs longer than {max_aa} aa ...", flush=True)
    with src.open() as inp, dest.open("w") as out:
        for line in inp:
            if line.startswith(">"):
                m = ORF_RE.match(line)
                length = int(m.group(5)) if m else max_aa + 1
                keep = length <= max_aa
                if keep:
                    kept += 1
                    out.write(line)
                else:
                    dropped += 1
            elif keep:
                out.write(line)
    print(f"Kept {kept} ORFs, dropped {dropped} giant stop-free stretches", flush=True)


def index_orfs(orf_fa: Path) -> dict[str, dict]:
    idx = {}
    with orf_fa.open() as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            m = ORF_RE.match(line)
            if not m:
                continue
            oid, src, c1, c2, length, frame = m.groups()
            idx[oid] = {
                "contig": src,
                "c1": int(c1),
                "c2": int(c2),
                "length": int(length),
                "frame": int(frame),
                "strand": "+" if int(c1) < int(c2) or int(frame) <= 3 else "-",
            }
            if idx[oid]["c1"] > idx[oid]["c2"]:
                idx[oid]["strand"] = "-"
            else:
                idx[oid]["strand"] = "+"
    print(f"Indexed {len(idx)} ORFs", flush=True)
    return idx


def aa_to_genome(orf: dict, ali_from: int, ali_to: int) -> tuple[int, int]:
    c1, c2 = orf["c1"], orf["c2"]
    if orf["strand"] == "+":
        g1 = c1 + (ali_from - 1) * 3
        g2 = c1 + ali_to * 3 - 1
        return min(g1, g2), max(g1, g2)
    g_hi = c1 - (ali_from - 1) * 3
    g_lo = c1 - (ali_to * 3 - 1)
    return min(g_lo, g_hi), max(g_lo, g_hi)


def genomic_hits(domtbl: Path, idx: dict, kind: str) -> list[dict]:
    hits = []
    if not domtbl.exists():
        return hits
    with domtbl.open() as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            col = line.split()
            score = float(col[13])
            if score < SCORE_MIN:
                continue
            oid = col[0]
            if oid not in idx:
                continue
            orf = idx[oid]
            ali_from, ali_to = int(col[17]), int(col[18])
            start, end = aa_to_genome(orf, ali_from, ali_to)
            hits.append(
                {
                    "kind": kind,
                    "contig": orf["contig"],
                    "strand": orf["strand"],
                    "start": start,
                    "end": end,
                    "score": score,
                    "orf": oid,
                    "ali_from": ali_from,
                    "ali_to": ali_to,
                }
            )
    return hits


def resolve_genomic(v_hits: list[dict], i_hits: list[dict]) -> list[dict]:
    cands = v_hits + i_hits
    cands.sort(key=lambda x: -x["score"])
    kept: list[dict] = []
    for c in cands:
        clash = False
        for k in kept:
            if c["contig"] != k["contig"] or c["strand"] != k["strand"]:
                continue
            ol = p.overlap_len(c["start"], c["end"], k["start"], k["end"])
            if ol == 0:
                continue
            shorter = min(c["end"] - c["start"] + 1, k["end"] - k["start"] + 1)
            if ol >= 120 or ol >= 0.5 * shorter:
                clash = True
                break
        if not clash:
            kept.append(c)
    return kept


def build_loci(hits: list[dict]) -> list[dict]:
    by = defaultdict(list)
    for h in hits:
        by[(h["contig"], h["strand"])].append(h)
    loci = []
    for (contig, strand), grp in by.items():
        grp.sort(key=lambda x: x["start"] if strand == "+" else -x["end"])
        used = [False] * len(grp)
        for i, h in enumerate(grp):
            if used[i] or h["kind"] != "V":
                continue
            used[i] = True
            loc = {
                "contig": contig,
                "strand": strand,
                "v": h,
                "i": None,
                "start": h["start"],
                "end": h["end"],
            }
            for j, o in enumerate(grp):
                if used[j] or o["kind"] != "I":
                    continue
                if strand == "+":
                    gap = o["start"] - h["end"]
                    order_ok = o["start"] >= h["start"]
                else:
                    gap = h["start"] - o["end"]
                    order_ok = o["end"] <= h["end"]
                if order_ok and 0 <= gap <= MAX_VI_GAP:
                    used[j] = True
                    loc["i"] = o
                    loc["start"] = min(loc["start"], o["start"])
                    loc["end"] = max(loc["end"], o["end"])
                    break
            loci.append(loc)
        for i, h in enumerate(grp):
            if used[i] or h["kind"] != "I":
                continue
            used[i] = True
            loci.append(
                {
                    "contig": contig,
                    "strand": strand,
                    "v": None,
                    "i": h,
                    "start": h["start"],
                    "end": h["end"],
                }
            )
    loci.sort(key=lambda x: (x["contig"], x["start"]))
    return loci


def write_windows(loci: list[dict], genome: dict[str, str], fa: Path) -> list[dict]:
    windows = []
    recs = []
    for i, loc in enumerate(loci, 1):
        ctg = loc["contig"]
        if ctg not in genome:
            continue
        clen = len(genome[ctg])
        w1 = max(1, loc["start"] - FLANK)
        w2 = min(clen, loc["end"] + FLANK)
        if w2 - w1 + 1 > SPAN_CAP:
            mid = (loc["start"] + loc["end"]) // 2
            w1 = max(1, mid - SPAN_CAP // 2)
            w2 = min(clen, w1 + SPAN_CAP - 1)
        wid = f"win_{i:04d}"
        recs.append((wid, genome[ctg][w1 - 1 : w2]))
        windows.append({**loc, "win_id": wid, "win_start": w1, "win_end": w2})
    p.write_fasta(fa, recs)
    print(f"Wrote {len(recs)} miniprot windows", flush=True)
    return windows


def map_window_models(gff: Path, windows: list[dict], genome: dict[str, str]) -> list[dict]:
    by_win = {w["win_id"]: w for w in windows}
    raw = p.parse_miniprot_gff(gff)
    models = []
    for i, g in enumerate(raw, 1):
        w = by_win.get(g["contig"])
        if w is None:
            continue
        off = w["win_start"] - 1
        cds = [{"start": e["start"] + off, "end": e["end"] + off, "phase": e["phase"]} for e in g["cds"]]
        gene = {
            "contig": w["contig"],
            "strand": g["strand"],
            "cds": cds,
        }
        nt, exons, notes = p.spliced_cds(genome, gene)
        aa = p.translate(nt) if nt else ""
        stops_mid = "*" in aa[:-1] if aa else True
        start = min(e["start"] for e in exons)
        end = max(e["end"] for e in exons)
        models.append(
            {
                "gene_id": f"hmm_{i:04d}",
                "mp_id": g["id"],
                "target": g["target"],
                "query": p.query_name(g["target"]),
                "contig": w["contig"],
                "start": start,
                "end": end,
                "strand": g["strand"],
                "span": end - start + 1,
                "n_exons": len(exons),
                "aa": aa,
                "nt_len": len(nt),
                "splice_ok": bool(nt) and not any(n.startswith("noncanonical") for n in notes),
                "orf_ok": bool(aa) and (not stops_mid) and len(aa) >= 80,
                "notes": notes,
                "exons": exons,
                "identity": g["identity"],
                "source": "hmm_window",
                "win_id": w["win_id"],
            }
        )
    return models


