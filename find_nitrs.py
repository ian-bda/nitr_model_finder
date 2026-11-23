#!/usr/bin/env python3
"""Find NITR loci in one reference genome.

Same search used for rainbow trout, Atlantic salmon, and northern pike.
Pass any nucleotide FASTA. Seeds default to the published proteins in seeds/.

Example:
  python3 find_nitrs.py --genome genome.fna.gz --name my_species --outdir results/my_species
"""
from __future__ import annotations

import argparse
import gzip
import os
import re
import shutil
import subprocess
import sys
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "nitrfind"))

import domain_annotation as ann  # noqa: E402
import domain_scan as h  # noqa: E402
import exon_repair as fix  # noqa: E402
import gene_models as p  # noqa: E402

DEFAULT_FULL = HERE / "seeds" / "round1_full.faa"
DEFAULT_V = HERE / "seeds" / "round1_V.faa"
DEFAULT_I = HERE / "seeds" / "round1_I.faa"

MINIPROT_OUTS = "0.25"
MINIPROT_OUTC = "0.15"
MINIPROT_G = "12k"
MAX_ORF_AA = 2000
OFFCLUSTER_V_MIN = 80.0
HYDRO = set("AILMFWV")
ITIM = re.compile(r"[IVL].Y..[LV]")
TOOLS = (
    "miniprot",
    "mafft",
    "hmmbuild",
    "hmmsearch",
    "hmmalign",
    "esl-translate",
    "samtools",
    "makeblastdb",
    "tblastn",
)


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Search one genome for NITR loci.")
    ap.add_argument("--genome", required=True, type=Path, help="Reference FASTA, plain or gzipped")
    ap.add_argument("--name", required=True, help="Short name used in locus ids and the report")
    ap.add_argument("--outdir", required=True, type=Path, help="Directory for models, tables, and proteins")
    ap.add_argument("--seeds", type=Path, default=DEFAULT_FULL, help="Full-length NITR proteins (FASTA)")
    ap.add_argument("--seeds-v", type=Path, default=DEFAULT_V, help="Extra V-domain peptides")
    ap.add_argument("--seeds-i", type=Path, default=DEFAULT_I, help="I-domain peptides. The I profile is built only from this file.")
    ap.add_argument(
        "--panel",
        type=Path,
        default=None,
        help="Optional full-length proteins from the target species or a close relative. Used as extra miniprot queries.",
    )
    ap.add_argument(
        "--check",
        type=Path,
        default=None,
        help="Optional held-out proteins. Compared with blastp after the search. Never used as queries.",
    )
    ap.add_argument("--threads", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", "8")))
    ap.add_argument("--rebuild-rounds", type=int, default=1, help="How many times to rebuild the V profile from accepted genes")
    return ap.parse_args()


def require_tools() -> None:
    missing = [name for name in TOOLS if shutil.which(name) is None]
    if missing:
        raise SystemExit(
            "Missing programs on PATH: "
            + ", ".join(missing)
            + "\nInstall with: conda install -c bioconda miniprot mafft hmmer samtools blast"
        )


def read_records(path: Path, min_len: int, max_len: int) -> list[tuple[str, str]]:
    if not path.exists():
        raise SystemExit(f"Missing FASTA: {path}")
    out = []
    seen = set()
    for hdr, aa in p.read_fasta(path).items():
        seq = aa.replace("*", "")
        if not (min_len <= len(seq) <= max_len) or seq in seen:
            continue
        seen.add(seq)
        out.append((p.fasta_id(hdr), seq))
    return out


def collect_seeds(args: argparse.Namespace) -> tuple[list[tuple[str, str]], list[tuple[str, str]], list[tuple[str, str]]]:
    full = read_records(args.seeds, 80, 10_000)
    n_seed = len(full)
    seen = {aa for _name, aa in full}
    if args.panel is not None:
        for name, aa in read_records(args.panel, 80, 10_000):
            if aa not in seen:
                seen.add(aa)
                full.append((name, aa))
    v_extra = read_records(args.seeds_v, 60, 160) if args.seeds_v.exists() else []
    i_dom = read_records(args.seeds_i, 60, 160)
    if not full:
        raise SystemExit(f"No full-length proteins in {args.seeds}")
    if not i_dom:
        raise SystemExit(f"No I-domain peptides in {args.seeds_i}")
    print(
        f"Seeds: {len(full)} full ({n_seed} from {args.seeds.name}"
        f"{'' if args.panel is None else f' + {len(full) - n_seed} panel'}), "
        f"{len(v_extra)} extra V, {len(i_dom)} I",
        flush=True,
    )
    return full, v_extra, i_dom


def prepare_genome(src: Path, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if not src.exists():
        raise SystemExit(f"Genome not found: {src}")
    if dest.exists() and dest.stat().st_size > 100_000_000:
        pass
    elif str(src).endswith(".gz"):
        print(f"Decompressing {src.name} ...", flush=True)
        with gzip.open(src, "rb") as inp, dest.open("wb") as out:
            while True:
                chunk = inp.read(8 * 1024 * 1024)
                if not chunk:
                    break
                out.write(chunk)
    else:
        if src.resolve() != dest.resolve():
            shutil.copyfile(src, dest)
    fai = Path(str(dest) + ".fai")
    if not fai.exists():
        p.run(["samtools", "faidx", str(dest)])
    return dest


def miniprot_gff(ref: Path, query: Path, dest: Path, threads: str) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"+ miniprot --outs {MINIPROT_OUTS} {ref.name}", flush=True)
    with dest.open("w") as fh:
        subprocess.check_call(
            [
                "miniprot",
                "-t",
                threads,
                "-j",
                "2",
                "-G",
                MINIPROT_G,
                "-N",
                "50",
                "--outs",
                MINIPROT_OUTS,
                "--outc",
                MINIPROT_OUTC,
                "--gff",
                str(ref),
                str(query),
            ],
            stdout=fh,
        )


def hmmsearch_pep(pep: Path, hmm_v: Path, hmm_i: Path, tag: str, out: Path, threads: str) -> tuple[dict, dict]:
    for kind, hmm in (("V", hmm_v), ("I", hmm_i)):
        tbl = out / "hmms" / f"{tag}_vs_NITR_{kind}.domtbl"
        aln = out / "hmms" / f"{tag}_vs_NITR_{kind}.out"
        print(f"+ hmmsearch NITR_{kind} ({tag})", flush=True)
        with aln.open("w") as fh:
            subprocess.check_call(
                ["hmmsearch", "--domtblout", str(tbl), "-E", "1e-5", "--cpu", threads, str(hmm), str(pep)],
                stdout=fh,
            )
    return (
        p.parse_hmmsearch_domtbl(out / "hmms" / f"{tag}_vs_NITR_V.domtbl"),
        p.parse_hmmsearch_domtbl(out / "hmms" / f"{tag}_vs_NITR_I.domtbl"),
    )


def models_from_gff(gff: Path, genome: dict[str, str], prefix: str, source: str) -> list[dict]:
    raw = p.parse_miniprot_gff(gff)
    models = []
    for i, g in enumerate(raw, 1):
        nt, exons, notes = p.spliced_cds(genome, g)
        aa = p.translate(nt) if nt else ""
        stops_mid = "*" in aa[:-1] if aa else True
        start = min(e["start"] for e in exons) if exons else g["start"]
        end = max(e["end"] for e in exons) if exons else g["end"]
        models.append(
            {
                "gene_id": f"{prefix}_{i:04d}",
                "mp_id": g["id"],
                "target": g["target"],
                "query": p.query_name(g["target"]),
                "contig": g["contig"],
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
                "source": source,
            }
        )
    return models


def annotate_protein(aa: str) -> dict:
    seq = (aa or "").replace("*", "")
    hits = list(p.G16C23.finditer(seq))
    tm = False
    tm_pos = 0
    for i in range(0, max(0, len(seq) - 19)):
        window = seq[i : i + 20]
        if sum(ch in HYDRO for ch in window) >= 15:
            tm = True
            tm_pos = i + 1
            break
    tail = seq[tm_pos + 19 :] if tm else seq[max(0, len(seq) - 80) :]
    return {
        "n_g16c23": len(hits),
        "g16_pos": ",".join(str(m.start() + 1) for m in hits) or ".",
        "tm": tm,
        "n_itim": len(ITIM.findall(tail)),
        "cys_n": seq.count("C"),
    }


def assign_status(row: dict) -> str:
    notes = ";".join(row.get("notes") or [])
    n_dom = int(row.get("n_v") or 0) + int(row.get("n_i") or 0)
    span = int(row.get("span") or 0)
    aa = (row.get("aa") or "").replace("*", "")
    if row.get("final_ok"):
        if n_dom == 1 and span < 1200:
            return "relic"
        if not aa.startswith("M"):
            return "truncated"
        return "complete"
    if "noncanonical" in notes:
        return "needs_review"
    if not row.get("orf_ok"):
        return "truncated"
    if n_dom == 0:
        return "truncated"
    if n_dom == 1 and span < 1500:
        return "relic"
    return "needs_review"


def qc_and_split(models, genome, hmm_v, hmm_i, tag, out, threads) -> list[dict]:
    if not models:
        return []
    pep = out / "proteins" / f"{tag}.raw.faa"
    p.write_fasta(pep, [(m["gene_id"], (m.get("aa") or "").replace("*", "")) for m in models])
    v_raw, i_raw = hmmsearch_pep(pep, hmm_v, hmm_i, tag + "_raw", out, threads)
    next_id = [len(models)]
    split_models = []
    for m in models:
        domains = p.resolve_domains(v_raw.get(m["gene_id"], []), i_raw.get(m["gene_id"], []))
        n_v = sum(1 for d in domains if d["kind"] == "V")
        n_i = sum(1 for d in domains if d["kind"] == "I")
        if n_v > 1 or n_i > 1:
            parts = p.split_fused_model(m, domains, genome, next_id)
            if len(parts) > 1:
                for part in parts:
                    part["source"] = m.get("source", tag)
                split_models.extend(parts)
                continue
        split_models.append(m)
    pep2 = out / "proteins" / f"{tag}.faa"
    p.write_fasta(pep2, [(m["gene_id"], (m.get("aa") or "").replace("*", "")) for m in split_models])
    v2, i2 = hmmsearch_pep(pep2, hmm_v, hmm_i, tag + "_split", out, threads)
    rows = []
    for m in split_models:
        domains = p.resolve_domains(v2.get(m["gene_id"], []), i2.get(m["gene_id"], []))
        n_v = sum(1 for d in domains if d["kind"] == "V")
        n_i = sum(1 for d in domains if d["kind"] == "I")
        fused = n_v > 1 or n_i > 1
        g16 = bool(p.G16C23.search(m.get("aa") or ""))
        span = m["end"] - m["start"] + 1
        if fused:
            m["notes"] = list(m.get("notes") or []) + [f"FUSED_V{n_v}_I{n_i}"]
        final_ok = (
            m.get("splice_ok", False)
            and m.get("orf_ok", False)
            and g16
            and not fused
            and span <= p.SPAN_CAP
            and n_v + n_i >= 1
        )
        v_hits = v2.get(m["gene_id"], [])
        i_hits = i2.get(m["gene_id"], [])
        v_best = max(v_hits, key=lambda x: x["score"], default=None)
        i_best = max(i_hits, key=lambda x: x["score"], default=None)

        def covers(hit: dict | None) -> bool:
            return bool(hit) and hit["hmm_from"] <= 16 and hit["hmm_to"] >= 23

        row = {
            **m,
            "n_v": n_v,
            "n_i": n_i,
            "fused": fused,
            "g16c23": g16,
            "span": span,
            "final_ok": final_ok,
            "v_score": v_best["score"] if v_best else 0,
            "i_score": i_best["score"] if i_best else 0,
            "v_aa_from": v_best["ali_from"] if v_best else "",
            "v_aa_to": v_best["ali_to"] if v_best else "",
            "i_aa_from": i_best["ali_from"] if i_best else "",
            "i_aa_to": i_best["ali_to"] if i_best else "",
            "v_g16c23_hmm": covers(v_best),
            "i_g16c23_hmm": covers(i_best),
        }
        row.update(annotate_protein(m.get("aa") or ""))
        ann.apply_leader_trim(row)
        row.update(annotate_protein(row.get("aa") or ""))
        row["status"] = assign_status(row)
        rows.append(row)
    return rows


TABLE_COLS = [
    "gene_id", "source", "status", "contig", "cluster", "start", "end", "strand", "span",
    "n_exons", "n_v", "n_i", "fused", "splice_ok", "orf_ok", "g16c23", "n_g16c23", "g16_pos",
    "v_g16c23_hmm", "i_g16c23_hmm", "tm", "n_itim", "cys_n", "final_ok", "query", "identity",
    "v_score", "i_score", "notes",
]


def write_table(path: Path, rows: list[dict], cluster_fn) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        fh.write("\t".join(TABLE_COLS) + "\n")
        for r in rows:
            notes = r.get("notes") or []
            note_s = notes if isinstance(notes, str) else (";".join(notes) or ".")
            vals = [
                r.get("gene_id", "."), r.get("source", "."), r.get("status", "."), r["contig"],
                cluster_fn(r["contig"]) if callable(cluster_fn) else r.get("cluster", "."),
                r["start"], r["end"], r.get("strand", "."), r.get("span", r["end"] - r["start"] + 1),
                r.get("n_exons", "."), r.get("n_v", "."), r.get("n_i", "."), r.get("fused", "."),
                r.get("splice_ok", "."), r.get("orf_ok", "."), r.get("g16c23", "."), r.get("n_g16c23", "."),
                r.get("g16_pos", "."), r.get("v_g16c23_hmm", "."), r.get("i_g16c23_hmm", "."),
                r.get("tm", "."), r.get("n_itim", "."), r.get("cys_n", "."), r.get("final_ok", "."),
                r.get("query", "."), r.get("identity", "."), r.get("v_score", "."), r.get("i_score", "."),
                note_s,
            ]
            fh.write("\t".join(str(v) for v in vals) + "\n")


def cluster_map(hmm_loci: list[dict]):
    counts = Counter(x["contig"] for x in hmm_loci if x.get("v") and x["v"]["score"] >= OFFCLUSTER_V_MIN)
    labels = {}
    ranked = 0
    for contig, n in counts.most_common():
        if n >= 3:
            ranked += 1
            labels[contig] = "cluster" if ranked == 1 else f"cluster{ranked}"
        elif n >= 1:
            labels[contig] = "hit"

    def lab(contig: str) -> str:
        return labels.get(contig, "other")

    return lab


def decorate_rows(rows: list[dict], lab) -> None:
    for r in rows:
        r["cluster"] = lab(r["contig"])
        r["status"] = assign_status(r)


def orf_search(genome_fa: Path, hmm_v: Path, hmm_i: Path, out: Path, threads: str) -> list[dict]:
    raw = out / "orfs" / "orfs.faa"
    filt = out / "orfs" / "orfs_max2k.faa"
    raw.parent.mkdir(exist_ok=True)
    if not raw.exists() or raw.stat().st_size < 1_000_000:
        print("esl-translate -l 60 -W ...", flush=True)
        with raw.open("w") as fh:
            subprocess.check_call(["esl-translate", "-l", "60", "-W", str(genome_fa)], stdout=fh)
    if not filt.exists() or filt.stat().st_size < 1_000_000:
        h.filter_orfs(raw, filt, max_aa=MAX_ORF_AA)
    idx = h.index_orfs(filt)
    for kind, hmm in (("V", hmm_v), ("I", hmm_i)):
        tbl = out / "hmms" / f"orfs_vs_NITR_{kind}.domtbl"
        aln = out / "hmms" / f"orfs_vs_NITR_{kind}.out"
        print(f"+ hmmsearch NITR_{kind} vs ORFs", flush=True)
        with aln.open("w") as fh:
            subprocess.check_call(
                ["hmmsearch", "--domtblout", str(tbl), "-E", "1e-5", "--cpu", threads, str(hmm), str(filt)],
                stdout=fh,
            )
    v_hits = h.genomic_hits(out / "hmms" / "orfs_vs_NITR_V.domtbl", idx, "V")
    i_hits = h.genomic_hits(out / "hmms" / "orfs_vs_NITR_I.domtbl", idx, "I")
    loci = h.build_loci(h.resolve_genomic(v_hits, i_hits))
    print(f"HMM loci: {len(loci)} (V hits {len(v_hits)}, I hits {len(i_hits)})", flush=True)
    return loci


def follow_locus(loc: dict) -> bool:
    vscore = loc["v"]["score"] if loc.get("v") else 0
    return vscore >= OFFCLUSTER_V_MIN or bool(loc.get("i"))


def uncovered_fn(complete: list[dict]):
    def uncovered(loc: dict) -> bool:
        for r in complete:
            if r["contig"] == loc["contig"] and p.overlap(loc["start"], loc["end"], r["start"], r["end"]):
                return False
        return True

    return uncovered


def held_out_check(unique: list[dict], out: Path, check: Path | None) -> None:
    dest = out / "validation"
    dest.mkdir(parents=True, exist_ok=True)
    if check is None:
        return
    faa = out / "proteins" / "unique_complete.faa"
    report = dest / "heldout_recovery.tsv"
    if not faa.exists() or faa.stat().st_size < 10:
        report.write_text("query\tstatus\n.\tno_catalog_proteins\n")
        return
    blast = dest / "heldout_vs_catalog.blastp"
    try:
        subprocess.check_call(
            ["makeblastdb", "-in", str(faa), "-dbtype", "prot", "-out", str(dest / "catalog")],
            stdout=subprocess.DEVNULL,
        )
        with blast.open("w") as fh:
            subprocess.check_call(
                [
                    "blastp", "-query", str(check), "-db", str(dest / "catalog"),
                    "-outfmt", "6 qseqid sseqid pident length qlen slen evalue bitscore",
                    "-evalue", "1e-5", "-max_target_seqs", "3",
                ],
                stdout=fh,
            )
    except (subprocess.CalledProcessError, FileNotFoundError) as exc:
        print(f"blastp check skipped: {exc}", flush=True)
        report.write_text("query\tstatus\n.\tblastp_unavailable\n")
        return
    best: dict[str, list[str]] = {}
    for line in blast.read_text().splitlines():
        col = line.split("\t")
        if len(col) < 8 or col[0] in best:
            continue
        ident = float(col[2])
        status = "found_complete" if ident >= 90 else "found_related" if ident >= 50 else "weak"
        best[col[0]] = [col[0], col[1], col[2], col[6], status]
    with report.open("w") as fh:
        fh.write("query\thit\tpident\tevalue\tstatus\n")
        for q in [p.fasta_id(h) for h in p.read_fasta(check)]:
            fh.write(("\t".join(best[q]) if q in best else f"{q}\t.\t.\t.\tnot_found") + "\n")
    print(f"Wrote {report}", flush=True)


def main() -> int:
    args = parse_args()
    require_tools()
    os.environ["SLURM_CPUS_PER_TASK"] = str(args.threads)
    p.THREADS = str(args.threads)
    h.THREADS = str(args.threads)
    threads = str(args.threads)
    gid = re.sub(r"[^A-Za-z0-9._-]", "_", args.name)
    out = args.outdir
    for d in ("genome", "seeds", "hmms", "miniprot", "proteins", "orfs", "tables", "validation"):
        (out / d).mkdir(parents=True, exist_ok=True)

    print(f"=== NITR search {gid} ===", flush=True)
    full, v_extra, i_dom = collect_seeds(args)
    full_faa, hmm_v, hmm_i = _build_hmms(full, v_extra, i_dom, out)
    genome_fa = prepare_genome(args.genome, out / "genome" / f"{gid}.fna")
    genome = p.load_genome(genome_fa)

    gff1 = out / "miniprot" / "seeds_vs_genome.gff"
    miniprot_gff(genome_fa, full_faa, gff1, threads)
    rows = qc_and_split(models_from_gff(gff1, genome, f"{gid}_mp", "miniprot"), genome, hmm_v, hmm_i, f"{gid}_mp", out, threads)

    loci = orf_search(genome_fa, hmm_v, hmm_i, out, threads)
    lab = cluster_map(loci)
    complete = [r for r in rows if r["final_ok"]]
    uncovered = uncovered_fn(complete)
    new_loci = [loc for loc in loci if follow_locus(loc) and uncovered(loc)]
    print(f"HMM loci not covered by miniprot complete models: {len(new_loci)}", flush=True)
    windows = h.write_windows(new_loci, genome, out / "miniprot" / "hmm_windows.fna")
    if windows:
        gff2 = out / "miniprot" / "seeds_vs_windows.gff"
        miniprot_gff(out / "miniprot" / "hmm_windows.fna", full_faa, gff2, threads)
        wmodels = h.map_window_models(gff2, windows, genome)
        for m in wmodels:
            m["gene_id"] = m["gene_id"].replace("hmm_", f"{gid}_w_")
            m["source"] = "hmm_window"
        rows.extend(qc_and_split(wmodels, genome, hmm_v, hmm_i, f"{gid}_w", out, threads))

    for rnd in range(args.rebuild_rounds):
        complete_aa, seen = [], set()
        for r in rows:
            if not r.get("final_ok"):
                continue
            aa = (r.get("aa") or "").replace("*", "")
            if len(aa) < 80 or aa in seen:
                continue
            seen.add(aa)
            complete_aa.append((r["gene_id"], aa))
        if len(complete_aa) < 8:
            print("Too few complete proteins to rebuild the V profile", flush=True)
            break
        print(f"Rebuild V profile round {rnd + 1} from {len(complete_aa)} accepted proteins + seeds", flush=True)
        accepted_v = []
        for name, aa in complete_aa:
            v, _i = p.split_v_i(aa)
            if v and 60 <= len(v) <= 160:
                accepted_v.append((name + "_V", v))
        full2, hmm_v, hmm_i = _build_hmms(full + complete_aa, v_extra + accepted_v, i_dom, out)
        filt = out / "orfs" / "orfs_max2k.faa"
        for kind, hmm in (("V", hmm_v), ("I", hmm_i)):
            tbl = out / "hmms" / f"orfs_rebuild_vs_NITR_{kind}.domtbl"
            aln = out / "hmms" / f"orfs_rebuild_vs_NITR_{kind}.out"
            print(f"+ hmmsearch NITR_{kind} vs ORFs (rebuild)", flush=True)
            with aln.open("w") as fh:
                subprocess.check_call(
                    ["hmmsearch", "--domtblout", str(tbl), "-E", "1e-5", "--cpu", threads, str(hmm), str(filt)],
                    stdout=fh,
                )
        idx = h.index_orfs(filt)
        v_hits = h.genomic_hits(out / "hmms" / "orfs_rebuild_vs_NITR_V.domtbl", idx, "V")
        i_hits = h.genomic_hits(out / "hmms" / "orfs_rebuild_vs_NITR_I.domtbl", idx, "I")
        extra = [
            loc
            for loc in h.build_loci(h.resolve_genomic(v_hits, i_hits))
            if follow_locus(loc) and uncovered_fn([r for r in rows if r.get("final_ok")])(loc)
        ]
        print(f"Rebuild extra loci: {len(extra)}", flush=True)
        if not extra:
            break
        windows2 = h.write_windows(extra, genome, out / "miniprot" / "rebuild_windows.fna")
        gff3 = out / "miniprot" / "rebuild_vs_windows.gff"
        miniprot_gff(out / "miniprot" / "rebuild_windows.fna", full2, gff3, threads)
        rm = h.map_window_models(gff3, windows2, genome)
        for m in rm:
            m["gene_id"] = m["gene_id"].replace("hmm_", f"{gid}_r{rnd + 1}_")
            m["source"] = "hmm_rebuild"
        rows.extend(qc_and_split(rm, genome, hmm_v, hmm_i, f"{gid}_r{rnd + 1}", out, threads))

    decorate_rows(rows, lab)
    unique = p.collapse_unique(rows)
    decorate_rows(unique, lab)

    panel_recs: list[tuple[str, str]] = []
    seen_aa: set[str] = set()
    if args.panel is not None:
        for name, seq in read_records(args.panel, 80, 10_000):
            if seq not in seen_aa:
                seen_aa.add(seq)
                panel_recs.append((name, seq))
    for r in unique:
        if r.get("status") != "complete":
            continue
        seq = (r.get("aa") or "").replace("*", "")
        if len(seq) >= 80 and seq.startswith("M") and seq not in seen_aa:
            seen_aa.add(seq)
            panel_recs.append((r["gene_id"] + "_remap", seq))
    if panel_recs:
        panel_fa = out / "seeds" / "remap_panel.faa"
        p.write_fasta(panel_fa, panel_recs)
        print(f"Remap queries: {len(panel_recs)} (panel + complete models)", flush=True)
        gff_t = out / "miniprot" / "remap_vs_genome.gff"
        miniprot_gff(genome_fa, panel_fa, gff_t, threads)
        tmodels = models_from_gff(gff_t, genome, f"{gid}_tr", "remap_miniprot")
        rows.extend(qc_and_split(tmodels, genome, hmm_v, hmm_i, f"{gid}_tr", out, threads))
        decorate_rows(rows, lab)
        unique = p.collapse_unique(rows)
        decorate_rows(unique, lab)

    panel_path = args.panel if args.panel is not None else out / "seeds" / "no_panel.faa"
    unique = fix.repair_loci(unique, genome, genome_fa, args.seeds_i, panel_path, out, gid)
    unique = qc_and_split(unique, genome, hmm_v, hmm_i, f"{gid}_fix", out, threads)
    decorate_rows(unique, lab)
    unique = p.collapse_unique(unique)
    decorate_rows(unique, lab)

    write_table(out / "tables" / "gene_table.tsv", rows, lab)
    write_table(out / "tables" / "unique_loci.tsv", unique, lab)
    complete_u = [r for r in unique if r.get("status") == "complete"]
    p.write_fasta(
        out / "proteins" / "unique_complete.faa",
        [(r["gene_id"], (r.get("aa") or "").replace("*", "")) for r in complete_u if r.get("aa")],
    )
    ann.write_exports(unique, out, v_hmm=hmm_v, i_hmm=hmm_i)
    held_out_check(unique, out, args.check)

    by_st = Counter(r.get("status", ".") for r in unique)
    by_cl = Counter(r.get("cluster", ".") for r in complete_u)
    lines = [
        f"# NITR search {gid}",
        "",
        f"Genome: {args.genome}",
        f"Seeds: {args.seeds}",
        f"I domains: {args.seeds_i}",
        f"Panel: {args.panel if args.panel else 'none'}",
        (
            f"The V profile was rebuilt from accepted proteins plus the original seeds "
            f"({args.rebuild_rounds} round{'s' if args.rebuild_rounds != 1 else ''})."
            if args.rebuild_rounds
            else "The V profile is the seed profile only."
        ),
        "The I profile is built from the original I peptides only.",
        f"miniprot --outs {MINIPROT_OUTS} --outc {MINIPROT_OUTC}; ORF max {MAX_ORF_AA} aa.",
        "Locus ids are search ids. Catalog names are a later step.",
        "",
        f"- unique loci: {len(unique)}  status {dict(by_st)}",
        f"- complete: {len(complete_u)}  by cluster {dict(by_cl)}",
        "",
        f"Proteins: {out / 'proteins' / 'unique_complete.faa'}",
        f"Loci: {out / 'tables' / 'unique_loci.tsv'}",
        f"Exons: {out / 'tables' / 'exons.tsv'}",
        f"Domains: {out / 'tables' / 'domains.tsv'}",
        "",
    ]
    (out / "REPORT.md").write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines), flush=True)
    return 0


def _build_hmms(full, v_extra, i_dom, out: Path):
    """MAFFT + hmmbuild. The I profile uses i_dom only."""
    v_recs = list(v_extra)
    seen = {aa for _n, aa in v_recs}
    for name, aa in full:
        v, _i = p.split_v_i(aa)
        if v and 60 <= len(v) <= 160 and v not in seen:
            seen.add(v)
            v_recs.append((name + "_V", v))
    full_fa = out / "seeds" / "seeds_full.faa"
    p.write_fasta(full_fa, full)
    p.write_fasta(out / "seeds" / "seeds_V.faa", v_recs)
    p.write_fasta(out / "seeds" / "seeds_I.faa", i_dom)
    print(f"HMM build: {len(v_recs)} V, {len(i_dom)} I", flush=True)
    for kind, recs in (("V", v_recs), ("I", i_dom)):
        fa = out / "seeds" / f"seeds_{kind}.faa"
        aln = out / "hmms" / f"NITR_{kind}.aln"
        hmm = out / "hmms" / f"NITR_{kind}.hmm"
        aln.parent.mkdir(parents=True, exist_ok=True)
        with aln.open("w") as fh:
            subprocess.check_call(["mafft", "--auto", "--amino", "--quiet", str(fa)], stdout=fh)
        p.run(["hmmbuild", "--amino", str(hmm), str(aln)])
    return full_fa, out / "hmms" / "NITR_V.hmm", out / "hmms" / "NITR_I.hmm"


if __name__ == "__main__":
    raise SystemExit(main())
