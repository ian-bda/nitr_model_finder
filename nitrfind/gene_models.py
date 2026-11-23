#!/usr/bin/env python3
"""Gene-model helpers for a NITR search.

Read and write FASTA, parse a miniprot GFF, splice exons, translate the coding
sequence, split a fused model into one locus per V or I domain, and collapse
overlapping hits.
"""
from __future__ import annotations

import gzip
import os
import re
import subprocess
from collections import defaultdict
from pathlib import Path

THREADS = str(os.environ.get("SLURM_CPUS_PER_TASK", "8"))

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")
G16C23 = re.compile(r"G.{6}C")


def read_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, list[str]] = {}
    header = None
    rec: list[str] = []
    fh_cm = gzip.open(path, "rt") if str(path).endswith(".gz") else path.open("rt")
    with fh_cm as fh:
        for line in fh:
            if line.startswith(">"):
                if header is not None:
                    seqs[header] = rec
                header = line[1:].strip()
                rec = []
            else:
                rec.append(line.strip())
        if header is not None:
            seqs[header] = rec
    return {h: "".join(s).replace(" ", "").upper() for h, s in seqs.items()}


def fasta_id(header: str) -> str:
    return header.split()[0]


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as out:
        for h, s in records:
            if not s:
                continue
            out.write(f">{h}\n")
            for i in range(0, len(s), 60):
                out.write(s[i : i + 60] + "\n")


def split_v_i(seq: str) -> tuple[str | None, str | None]:
    """First GxxxxxxC after pos 10 is V; a second well-spaced one is I."""
    hits = [m.start() for m in G16C23.finditer(seq) if m.start() >= 10]
    v = i_dom = None
    if hits:
        v_start = max(0, hits[0] - 15)
        v_end = min(len(seq), v_start + 130)
        v = seq[v_start:v_end]
        later = [i for i in hits[1:] if i >= v_start + 80]
        if later:
            i_start = max(0, later[0] - 15)
            i_end = min(len(seq), i_start + 130)
            i_dom = seq[i_start:i_end]
    return v, i_dom


def run(cmd: list[str]) -> None:
    print("+", " ".join(cmd), flush=True)
    subprocess.check_call(cmd)


def load_genome(fa: Path) -> dict[str, str]:
    print("Loading genome into memory...", flush=True)
    recs = {}
    name = None
    chunks: list[str] = []
    with fa.open() as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    recs[name] = "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip().upper())
        if name is not None:
            recs[name] = "".join(chunks)
    print(f"Loaded {len(recs)} contigs", flush=True)
    return recs


def parse_miniprot_gff(path: Path) -> list[dict]:
    genes: dict[str, dict] = {}
    with path.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            seqid, _src, ftype, start, end, score, strand, phase, attr = line.rstrip("\n").split("\t")
            attrs = dict(x.split("=", 1) for x in attr.split(";") if "=" in x)
            if ftype == "mRNA":
                gid = attrs.get("ID", "")
                genes[gid] = {
                    "id": gid,
                    "contig": seqid,
                    "start": int(start),
                    "end": int(end),
                    "strand": strand,
                    "score": score,
                    "target": attrs.get("Target", ""),
                    "identity": attrs.get("Identity", attrs.get("identity", "")),
                    "cds": [],
                }
            elif ftype == "CDS":
                parent = attrs.get("Parent", "")
                if parent in genes:
                    genes[parent]["cds"].append(
                        {
                            "start": int(start),
                            "end": int(end),
                            "phase": phase if phase != "." else "0",
                        }
                    )
    out = []
    for g in genes.values():
        g["cds"].sort(key=lambda x: x["start"])
        if g["strand"] == "-":
            g["cds"] = list(reversed(g["cds"]))
        if g["cds"]:
            out.append(g)
    return out


def slice_seq(genome: dict[str, str], contig: str, start: int, end: int) -> str:
    seq = genome[contig]
    if start < 1 or end > len(seq) or end < start:
        return ""
    return seq[start - 1 : end]


def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def spliced_cds(genome: dict[str, str], gene: dict) -> tuple[str, list[dict], list[str]]:
    notes = []
    parts = []
    exons_out = []
    cds = gene["cds"]
    if gene["contig"] not in genome:
        return "", [], ["missing_contig"]
    for i, ex in enumerate(cds):
        seq = slice_seq(genome, gene["contig"], ex["start"], ex["end"])
        if gene["strand"] == "-":
            seq = revcomp(seq)
        parts.append(seq)
        donor = acc = "."
        if i < len(cds) - 1:
            nxt = cds[i + 1]
            if gene["strand"] == "+":
                donor = slice_seq(genome, gene["contig"], ex["end"] + 1, ex["end"] + 2)
                acc = slice_seq(genome, gene["contig"], nxt["start"] - 2, nxt["start"] - 1)
            else:
                donor = revcomp(slice_seq(genome, gene["contig"], ex["start"] - 2, ex["start"] - 1))
                acc = revcomp(slice_seq(genome, gene["contig"], nxt["end"] + 1, nxt["end"] + 2))
            if donor != "GT" or acc != "AG":
                notes.append(f"noncanonical_junc_{i}_{donor}..{acc}")
        exons_out.append({**ex, "seq_len": len(seq), "donor": donor, "acceptor": acc})
    return "".join(parts), exons_out, notes


def translate(nt: str) -> str:
    table = {
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
        "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
        "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
        "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    }
    aa = []
    for i in range(0, len(nt) - 2, 3):
        aa.append(table.get(nt[i : i + 3], "X"))
    return "".join(aa)


def parse_hmmsearch_domtbl(path: Path) -> dict[str, list[dict]]:
    hits: dict[str, list[dict]] = defaultdict(list)
    if not path.exists():
        return hits
    with path.open() as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            p = line.split()
            hits[p[0]].append(
                {
                    "query": p[3],
                    "evalue": float(p[12]),
                    "score": float(p[13]),
                    "hmm_from": int(p[15]),
                    "hmm_to": int(p[16]),
                    "ali_from": int(p[17]),
                    "ali_to": int(p[18]),
                }
            )
    return hits


SCORE_MIN = 20.0
SPAN_CAP = 15_000


def overlap(a1: int, a2: int, b1: int, b2: int) -> bool:
    return a1 <= b2 and b1 <= a2


def overlap_len(a1: int, a2: int, b1: int, b2: int) -> int:
    return max(0, min(a2, b2) - max(a1, b1) + 1)


def same_domain_region(a: dict, b: dict) -> bool:
    """True if two HMM hits are the same Ig exon, not V then I next door.

    A few aa of overlap at the V/I junction is normal and must not collapse
    the real V into the I call.
    """
    ol = overlap_len(a["ali_from"], a["ali_to"], b["ali_from"], b["ali_to"])
    if ol == 0:
        return False
    shorter = min(a["ali_to"] - a["ali_from"] + 1, b["ali_to"] - b["ali_from"] + 1)
    return ol >= 40 or ol >= 0.5 * shorter


def query_name(target: str) -> str:
    return target.split()[0] if target else "."


def resolve_domains(v_hits: list[dict], i_hits: list[dict]) -> list[dict]:
    """One call per region. If V and I overlap, keep the higher bit score.

    The V HMM often hits the I exon (and vice versa). That is not two genes.
    """
    cands = [{**h, "kind": "V"} for h in v_hits if h["score"] >= SCORE_MIN]
    cands += [{**h, "kind": "I"} for h in i_hits if h["score"] >= SCORE_MIN]
    cands.sort(key=lambda x: -x["score"])
    kept: list[dict] = []
    for c in cands:
        if any(same_domain_region(c, k) for k in kept):
            continue
        kept.append(c)
    kept.sort(key=lambda x: x["ali_from"])
    return kept


def split_exons_at_aa(exons: list[dict], strand: str, last_aa: int) -> tuple[list[dict], list[dict]]:
    """Split CDS exons after last_aa (1-based) of the spliced protein."""
    cut_nt = max(0, last_aa * 3)
    acc = 0
    left: list[dict] = []
    right: list[dict] = []
    for ex in exons:
        n = ex["seq_len"]
        if acc + n <= cut_nt:
            left.append(dict(ex))
        elif acc >= cut_nt:
            right.append(dict(ex))
        else:
            take = cut_nt - acc
            rest = n - take
            if take <= 0:
                right.append(dict(ex))
            elif rest <= 0:
                left.append(dict(ex))
            elif strand == "+":
                left_end = ex["start"] + take - 1
                left.append({**ex, "end": left_end, "seq_len": take, "donor": "."})
                right.append({**ex, "start": left_end + 1, "seq_len": rest, "acceptor": "."})
            else:
                left_start = ex["end"] - take + 1
                left.append({**ex, "start": left_start, "seq_len": take, "donor": "."})
                right.append({**ex, "end": left_start - 1, "seq_len": rest, "acceptor": "."})
        acc += n
    left = [e for e in left if e["end"] >= e["start"] and e["seq_len"] > 0]
    right = [e for e in right if e["end"] >= e["start"] and e["seq_len"] > 0]
    return left, right


def rebuild_from_exons(parent: dict, exons: list[dict], genome: dict[str, str], gene_id: str, extra_notes: list[str]) -> dict:
    fake = {"contig": parent["contig"], "strand": parent["strand"], "cds": exons}
    nt, exons2, splice_notes = spliced_cds(genome, fake)
    aa = translate(nt) if nt else ""
    stops_mid = "*" in aa[:-1] if aa else True
    start = min(e["start"] for e in exons2)
    end = max(e["end"] for e in exons2)
    return {
        "gene_id": gene_id,
        "mp_id": parent.get("mp_id", ""),
        "target": parent.get("target", ""),
        "query": parent.get("query", "."),
        "contig": parent["contig"],
        "start": start,
        "end": end,
        "strand": parent["strand"],
        "span": end - start + 1,
        "n_exons": len(exons2),
        "aa": aa,
        "nt_len": len(nt),
        "splice_ok": bool(nt) and not any(n.startswith("noncanonical") for n in splice_notes),
        "orf_ok": bool(aa) and (not stops_mid) and len(aa) >= 80,
        "notes": extra_notes + splice_notes,
        "exons": exons2,
        "identity": parent.get("identity", ""),
        "split_from": parent.get("gene_id", ""),
    }


def locus_cuts(domains: list[dict], aa_len: int) -> list[int] | None:
    """Return last-aa of each locus except the final one, or None if no split."""
    vs = [d for d in domains if d["kind"] == "V"]
    ihs = [d for d in domains if d["kind"] == "I"]
    anchors = vs if len(vs) > 1 else (ihs if len(ihs) > 1 else [])
    if len(anchors) < 2:
        return None
    cuts = []
    for i, a in enumerate(anchors[:-1]):
        nxt = anchors[i + 1]
        cut = nxt["ali_from"] - 1
        if a["ali_to"] >= nxt["ali_from"]:
            cut = (a["ali_to"] + nxt["ali_from"]) // 2
        cut = min(max(cut, 1), aa_len - 1)
        cuts.append(cut)
    return cuts


def split_fused_model(model: dict, domains: list[dict], genome: dict[str, str], next_id: list[int]) -> list[dict]:
    cuts = locus_cuts(domains, len(model["aa"]))
    if not cuts:
        return [model]
    remaining = list(model["exons"])
    consumed_aa = 0
    pieces: list[dict] = []
    for i, cut in enumerate(cuts, 1):
        local_cut = cut - consumed_aa
        if local_cut <= 0 or not remaining:
            break
        left, right = split_exons_at_aa(remaining, model["strand"], local_cut)
        if not left or not right:
            break
        next_id[0] += 1
        gid = f"pilot_{next_id[0]:04d}"
        part = rebuild_from_exons(
            model, left, genome, gid, [f"SPLIT_FROM_{model['gene_id']}_part{i}"]
        )
        pieces.append(part)
        consumed_aa += max(1, part["nt_len"] // 3)
        remaining = right
    if remaining:
        next_id[0] += 1
        gid = f"pilot_{next_id[0]:04d}"
        pieces.append(
            rebuild_from_exons(
                model,
                remaining,
                genome,
                gid,
                [f"SPLIT_FROM_{model['gene_id']}_part{len(pieces) + 1}"],
            )
        )
    return pieces if len(pieces) >= 2 else [model]


def _final_ok(row: dict) -> bool:
    v = row.get("final_ok")
    return v is True or v == "True"


def _identity(row: dict) -> float:
    try:
        return float(row.get("identity") or 0)
    except (TypeError, ValueError):
        return 0.0


def _aa_len(row: dict) -> int:
    aa = (row.get("aa") or "").replace("*", "")
    if aa:
        return len(aa)
    try:
        return int(row.get("nt_len") or 0) // 3
    except (TypeError, ValueError):
        return 0


def drop_spanning_artifacts(rows: list[dict]) -> tuple[list[dict], int]:
    """Drop a model that overlaps two non-overlapping complete models.

    miniprot often stitches tandem NITR paralogs into one long alignment.
    Those spans are not genes; keeping the longest one hides the real copies.
    """
    by = defaultdict(list)
    for i, r in enumerate(rows):
        by[(r["contig"], r["strand"])].append(i)
    drop: set[int] = set()
    for idxs in by.values():
        anchors = [i for i in idxs if _final_ok(rows[i]) and int(rows[i].get("span") or 0) >= 400]
        for i in idxs:
            b = rows[i]
            bspan = int(b.get("span") or 0)
            hits = [
                j
                for j in anchors
                if j != i
                and int(rows[j].get("span") or 0) < bspan
                and overlap(b["start"], b["end"], rows[j]["start"], rows[j]["end"])
            ]
            found = False
            for a, x in enumerate(hits):
                for y in hits[a + 1 :]:
                    if not overlap(rows[x]["start"], rows[x]["end"], rows[y]["start"], rows[y]["end"]):
                        drop.add(i)
                        found = True
                        break
                if found:
                    break
    kept = [r for i, r in enumerate(rows) if i not in drop]
    return kept, len(drop)


def collapse_unique(rows: list[dict]) -> list[dict]:
    """One row per non-overlapping locus on a contig+strand.

    Spanning tandem stitches are removed first. Among remaining overlaps
    (same gene, different seeds) keep complete, then more domains, more
    exons, longer protein. A short perfect fragment must not beat a full gene.
    """
    rows, n_drop = drop_spanning_artifacts(rows)
    if n_drop:
        print(
            f"Dropped {n_drop} spanning miniprot artifacts "
            "(overlap two non-overlapping complete models)",
            flush=True,
        )
    by = defaultdict(list)
    for r in rows:
        by[(r["contig"], r["strand"])].append(r)
    uniq = []
    for grp in by.values():
        grp = sorted(grp, key=lambda x: x["start"])
        used = [False] * len(grp)
        for i, a in enumerate(grp):
            if used[i]:
                continue
            cluster = [a]
            used[i] = True
            changed = True
            while changed:
                changed = False
                for j, b in enumerate(grp):
                    if used[j]:
                        continue
                    if any(overlap(b["start"], b["end"], c["start"], c["end"]) for c in cluster):
                        cluster.append(b)
                        used[j] = True
                        changed = True
            cluster.sort(
                key=lambda x: (
                    not _final_ok(x),
                    -(int(x.get("n_v") or 0) + int(x.get("n_i") or 0)),
                    -int(x.get("n_exons") or 0),
                    -_aa_len(x),
                    -_identity(x),
                    int(x.get("span") or 10**9),
                    x.get("gene_id", ""),
                )
            )
            rep = dict(cluster[0])
            rep["members"] = ",".join(x["gene_id"] for x in cluster)
            rep["n_members"] = len(cluster)
            syms = sorted(
                {
                    s
                    for x in cluster
                    for s in str(x.get("known") or ".").split(",")
                    if s not in {".", "", "None"}
                }
            )
            rep["known"] = ",".join(syms) or "."
            uniq.append(rep)
    uniq.sort(key=lambda x: (x["contig"], x["start"]))
    return uniq


