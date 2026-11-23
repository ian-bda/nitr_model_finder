#!/usr/bin/env python3
"""Place NITR domains and write the exon, domain, and GFF tables.

Each protein is labeled leader, V, I, transmembrane, and cytoplasmic tail.
"""
from __future__ import annotations

import re
import subprocess
from pathlib import Path

import gene_models as p

HYDRO = set("AILMFWV")
ITIM = re.compile(r"[IVL].Y..[LV]")
ITIM_TYPE = re.compile(r"[IVL].Y..[IVL]")
ITSM = re.compile(r"T.Y..[IV]")
CHARGED_TM = set("KR")
IG_MAX = 130
I_AFTER_V = 50
I_MIN_SCORE = 50.0
I_MIN_LEN = 70
LEADER_MIN = 10
LEADER_MAX = 40
I_MOTIFS = ("VTLNC", "VTLQC", "VALNC", "VSLNC", "LTLNC", "MTLQC", "VTLTC")
DISPLAY = {
    "leader": "SP",
    "V": "V",
    "I": "I",
    "TM": "TM",
    "CYT": "CT",
    "CYT1": "CT",
    "CYT2": "CT",
}


def find_tm(seq: str, after: int = 0) -> tuple[int, int] | None:
    """Best 20-aa hydrophobic window after `after` (0-based). 1-based output."""
    best = None
    best_n = 0
    start = max(0, after)
    for i in range(start, max(start, len(seq) - 19)):
        w = seq[i : i + 20]
        n = sum(ch in HYDRO for ch in w)
        if n >= 12 and n > best_n:
            best_n = n
            best = (i + 1, i + 20)
    return best


def _as_hits(x) -> list[dict]:
    if not x:
        return []
    return x if isinstance(x, list) else [x]


def _cap(start: int, end: int, max_len: int = IG_MAX) -> tuple[int, int]:
    if start <= 0 or end < start:
        return 0, 0
    return start, min(end, start + max_len - 1)


def pick_v_hit(v_hits, seq_len: int) -> dict | None:
    hits = _as_hits(v_hits)
    if not hits:
        return None
    half = max(40, seq_len // 2)
    early = [h for h in hits if int(h["ali_from"]) <= half]
    return max(early or hits, key=lambda h: float(h.get("score") or 0))


def pick_i_hit(i_hits, v_from: int, v_to: int, n_i, seq_len: int, tm=None) -> dict | None:
    """I must sit after V and before TM. Weak / TM-overlapping hits are not I."""
    hits = _as_hits(i_hits)
    after = (v_to + 1) if v_to else (v_from + I_AFTER_V if v_from else 0)
    tm_start = int(tm[0]) if tm else seq_len + 1
    kept = []
    for h in hits:
        if float(h.get("score") or 0) < I_MIN_SCORE:
            continue
        a = max(int(h["ali_from"]), after)
        b = min(int(h["ali_to"]), tm_start - 1, seq_len)
        if b - a + 1 < I_MIN_LEN:
            continue
        kept.append({**h, "ali_from": a, "ali_to": min(b, a + IG_MAX - 1)})
    if not kept:
        return None
    return max(kept, key=lambda h: float(h.get("score") or 0))


def ig_motif(seq: str, cys_pos: int) -> str:
    """Five residues ending at the cysteine."""
    if cys_pos < 5:
        return ""
    return seq[cys_pos - 5 : cys_pos]


def finish_ig_ends(seq: str, v: dict | None, i_hit: dict | None) -> None:
    if v:
        second = v.get("second") or 104
        v_end = min(len(seq), v["start"] + second + 11, v["start"] + 119)
        if i_hit:
            v_end = min(v_end, i_hit["start"] - 1)
        v["end"] = max(v_end, v["c23"])
    if i_hit:
        last = i_hit.get("last_cys") or 23
        i_end = min(len(seq), i_hit["start"] + last + 11, i_hit["start"] + 119)
        i_hit["end"] = max(i_end, i_hit["c23"])


def truncated_i(seq: str) -> dict | None:
    """An I that starts at residue 1 because the front of the domain is missing."""
    head = seq[:40]
    for motif in ("VTLNC", "VTLQC", "VTLTC", "PGDS"):
        at = head.find(motif)
        if at < 0:
            continue
        if motif == "PGDS":
            tail = head[at:]
            c_at = tail.find("C")
            if c_at < 0:
                continue
            c23 = at + c_at + 1
        else:
            c23 = at + len(motif)
        return {
            "start": 1,
            "end": min(len(seq), 115),
            "c23": c23,
            "g16": len(seq) >= 16 and seq[15] == "G",
            "n_cys": seq[: min(len(seq), 115)].count("C"),
            "last_cys": 23,
            "truncated": True,
            "source": "cys_frame",
        }
    return None


def place_ig_frames(seq: str) -> tuple[dict | None, dict | None]:
    """Place V and I so position 23 is cysteine.

    V: G16, no C30, a later Cys around the far end. A V domain often has
    T or L before C23, so that motif is not required.
    I: C30 plus extra Cys, or a PGDS / VTLNC motif in front of C23.
    """
    seq = (seq or "").replace("*", "")
    if not seq:
        return None, None
    cys = [i for i, aa in enumerate(seq, start=1) if aa == "C"]
    v_hits: list[dict] = []
    i_hits: list[dict] = []
    for cys_pos in cys:
        start = cys_pos - 22
        if start < 1:
            continue
        end_lim = min(len(seq), start + 114)
        motif = ig_motif(seq, cys_pos)
        window = seq[max(0, cys_pos - 30) : cys_pos]
        canonical_i = motif in I_MOTIFS or "PGDS" in window
        if end_lim - start + 1 < (40 if canonical_i else 60):
            continue
        g16 = start + 14 <= len(seq) and seq[start + 14] == "G"
        c30 = start + 29 <= len(seq) and seq[start + 28] == "C"
        n_cys = seq[start - 1 : end_lim].count("C")
        later = [pos - start + 1 for pos in cys if start < pos <= end_lim]
        frame = {
            "start": start,
            "c23": cys_pos,
            "g16": g16,
            "n_cys": n_cys,
            "last_cys": max(later) if later else 23,
            "canonical": canonical_i,
            "motif": motif,
            "source": "cys_frame",
        }
        if canonical_i and 2 <= n_cys <= 8:
            i_hits.append(frame)
        elif c30 and 4 <= n_cys <= 8:
            i_hits.append(frame)
        elif (
            g16
            and (not c30)
            and 1 <= n_cys <= 5
            and any(40 <= pos <= 120 for pos in later)
        ):
            frame["second"] = min(pos for pos in later if 40 <= pos <= 120)
            v_hits.append(frame)
    v_hits.sort(key=lambda hit: (not hit["g16"], hit["start"]))
    i_hits.sort(key=lambda hit: hit["start"])
    g16_hits = [hit for hit in v_hits if hit["g16"]]
    v = (g16_hits or v_hits or [None])[0]
    after = v["start"] + 50 if v else 0
    later_i = [hit for hit in i_hits if hit["start"] >= after]
    canonical = [hit for hit in later_i if hit.get("canonical")]
    i_hit = (canonical or later_i or [None])[0]
    if i_hit is None and v is None:
        i_hit = truncated_i(seq)
    finish_ig_ends(seq, v, i_hit if i_hit and not i_hit.get("truncated") else None)
    if i_hit and i_hit.get("truncated") and "end" not in i_hit:
        i_hit["end"] = min(len(seq), 115)
    return v, i_hit


def hmm_hit_to_frame(hit: dict | None, seq: str, kind: str) -> dict | None:
    if not hit:
        return None
    start = int(hit["ali_from"])
    end = int(hit["ali_to"])
    if end < start:
        return None
    c23 = start + 22 if start + 21 <= len(seq) else start
    if 1 <= c23 <= len(seq) and seq[c23 - 1] != "C":
        nearby = [
            i
            for i, aa in enumerate(seq[start - 1 : min(len(seq), start + 40)], start=start)
            if aa == "C"
        ]
        if nearby:
            c23 = min(nearby, key=lambda pos: abs(pos - (start + 22)))
    return {
        "start": start,
        "end": min(end, start + IG_MAX - 1, len(seq)),
        "c23": c23,
        "g16": start + 14 <= len(seq) and seq[start + 14] == "G",
        "n_cys": seq[start - 1 : end].count("C"),
        "last_cys": 23,
        "source": "hmm",
        "score": hit.get("score", "."),
        "kind": kind,
    }


def choose_ig(seq: str, v_hit: dict | None = None, i_hit: dict | None = None) -> tuple[dict | None, dict | None]:
    """Cysteine-frame V/I first; HMM fills a domain the frame missed."""
    seq = (seq or "").replace("*", "")
    v_cys, i_cys = place_ig_frames(seq)
    v_hmm = hmm_hit_to_frame(v_hit, seq, "V") if v_hit else None
    i_hmm = hmm_hit_to_frame(i_hit, seq, "I") if i_hit else None
    v = v_hmm or v_cys
    i = i_cys or i_hmm
    if v and i and i["start"] < v["start"] + 40:
        i = i_hmm if (i_hmm and i_hmm["start"] >= v["start"] + 40) else None
    if v and i and i["start"] <= v["end"]:
        v = dict(v)
        v["end"] = max(v["c23"], i["start"] - 1)
    return v, i


def receptor_tms(seq: str, v: dict | None, i: dict | None) -> list[tuple[int, int]]:
    """Hydrophobic spans after the last Ig domain. The Ig domain is trimmed back."""
    last = i or v
    if not last:
        return []
    ig_end = last["end"]
    tm = find_tm(seq, max(0, ig_end - 20))
    if not tm or tm[0] < last["start"] + 50:
        return []
    if tm[0] <= last["end"]:
        last["end"] = max(last["c23"], tm[0] - 1)
    spans = [tm]
    second = find_tm(seq, tm[1])
    if second and second[0] > tm[1]:
        spans.append(second)
    return spans


def signaling_type(seq: str, tms: list[tuple[int, int]], ig_end: int) -> str:
    """Inhibitory, activating, or both, from the ITIM and charged-transmembrane patterns."""
    if len(tms) > 1:
        tm_start, tm_end = tms[-1]
        tail = seq[tm_end:]
        tm_seq = seq[tm_start - 1 : tm_end]
        found = []
        if ITIM_TYPE.search(tail):
            found.append("ITIM")
        if ITSM.search(tail):
            found.append("ITSM")
        if any(residue in CHARGED_TM for residue in tm_seq):
            found.append("K/R in the last span")
        note = "; " + " and ".join(found) if found else ""
        return f"Ambiguous: more than one TM span{note}"
    if len(tms) == 1:
        tm_start, tm_end = tms[0]
        tail = seq[tm_end:]
        tm_seq = seq[tm_start - 1 : tm_end]
        charged = any(residue in CHARGED_TM for residue in tm_seq)
    elif ig_end:
        tail = seq[ig_end:]
        charged = False
    else:
        tail = ""
        charged = False
    has_itim = bool(ITIM_TYPE.search(tail))
    has_itsm = bool(ITSM.search(tail))
    found = []
    if has_itim:
        found.append("ITIM")
    if has_itsm:
        found.append("ITSM")
    if charged:
        found.append("K/R in TM")
    if has_itim and not charged:
        motifs = "ITIM and ITSM" if has_itsm else "ITIM"
        return f"Inhibiting: {motifs}"
    if charged and not has_itim and not has_itsm:
        return "Activating: K/R in TM"
    if has_itsm and not has_itim and not charged:
        return "Switch: ITSM"
    if found:
        return "Ambiguous: " + " and ".join(found)
    return "Ambiguous: no ITIM, ITSM, or charged TM"


def overlap_len(a1: int, a2: int, b1: int, b2: int) -> int:
    lo = max(a1, b1)
    hi = min(a2, b2)
    if hi < lo:
        return 0
    return hi - lo + 1


def exon_aa_bounds(exon: dict, seq_len: int) -> tuple[int, int] | None:
    try:
        start = int(exon["aa_start"])
        end = int(exon["aa_end"])
    except (TypeError, ValueError, KeyError):
        return None
    start = max(1, start)
    end = min(seq_len, end)
    if end < start:
        return None
    return start, end


def leader_call(exons: list[dict], seq_len: int, dom_start: int) -> tuple[bool, list[str], int]:
    """A candidate leader is one short exon that ends before the Ig domain."""
    flags: list[str] = []
    if dom_start <= 1:
        return False, flags, 0
    before = []
    holder = None
    for exon in exons:
        bounds = exon_aa_bounds(exon, seq_len)
        if not bounds:
            continue
        a, b = bounds
        if b < dom_start:
            before.append((b - a + 1, b))
        elif a <= dom_start <= b:
            holder = a
    if len(before) == 1 and LEADER_MIN <= before[0][0] <= LEADER_MAX:
        return True, flags, before[0][1]
    if len(before) > 1:
        flags.append("leader split across exons")
    elif len(before) == 1 and before[0][0] > LEADER_MAX:
        flags.append("long exon before the domain")
    elif holder and dom_start - holder >= LEADER_MIN:
        flags.append("leader shares the domain exon")
    return False, flags, 0


def regions_of(
    seq: str,
    v: dict | None,
    i: dict | None,
    tms: list[tuple[int, int]],
    leader_end: int,
) -> list[tuple[str, int, int]]:
    regions = []
    if leader_end >= 1:
        regions.append(("SP", 1, leader_end))
    if v:
        regions.append(("V", v["start"], v["end"]))
    if i:
        regions.append(("I", i["start"], i["end"]))
    if tms:
        regions.append(("TM", tms[0][0], tms[0][1]))
        if tms[0][1] < len(seq):
            regions.append(("CT", tms[0][1] + 1, len(seq)))
    return regions


def exon_call(seq: str, start: int, end: int, regions: list[tuple[str, int, int]]) -> tuple[str, str]:
    hits = []
    for name, r1, r2 in regions:
        n = overlap_len(start, end, r1, r2)
        if n:
            hits.append((name, r1, r2, n))
    use = [hit for hit in hits if hit[3] >= 10]
    if not use:
        return "", seq[start - 1 : end]
    if len(use) >= 2:
        label = "+".join(hit[0] for hit in use)
        lo = min(max(start, hit[1]) for hit in use)
        hi = max(min(end, hit[2]) for hit in use)
        return label, seq[lo - 1 : hi]
    name, r1, r2, _n = use[0]
    lo = max(start, r1)
    hi = min(end, r2)
    return name, seq[lo - 1 : hi]


def exon_at(exons: list[dict], seq_len: int, pos: int) -> int | None:
    for i, exon in enumerate(exons, start=1):
        bounds = exon_aa_bounds(exon, seq_len)
        if bounds and bounds[0] <= pos <= bounds[1]:
            try:
                return int(exon.get("exon") or i)
            except (TypeError, ValueError):
                return i
    return None


def exon_display(seq: str, exon: dict, exons: list[dict], v: dict | None, i: dict | None, tms: list, leader_end: int) -> tuple[str, str]:
    """V and I peptides start 22 residues before their cysteine.

    The peptide is written on the exon that contains that cysteine. A fused
    V+I+TM exon is labeled V+I, not TM.
    """
    bounds = exon_aa_bounds(exon, len(seq))
    if not bounds:
        return "", ""
    try:
        num = int(exon.get("exon") or 0)
    except (TypeError, ValueError):
        num = 0
    if not num:
        num = exons.index(exon) + 1 if exon in exons else 1
    framed = []
    for kind, dom in (("V", v), ("I", i)):
        if not dom:
            continue
        holder = exon_at(exons, len(seq), dom["c23"])
        if holder is None:
            holder = exon_at(exons, len(seq), dom["start"])
        if holder == num:
            framed.append((kind, seq[dom["start"] - 1 : dom["end"]]))
    if framed:
        return "+".join(kind for kind, _pep in framed), "\n".join(pep for _kind, pep in framed)
    regions = [
        region
        for region in regions_of(seq, v, i, tms, leader_end)
        if region[0] not in {"V", "I"}
    ]
    return exon_call(seq, bounds[0], bounds[1], regions)


def frame_to_hit(frame: dict | None) -> dict | None:
    if not frame:
        return None
    return {
        "ali_from": frame["start"],
        "ali_to": frame["end"],
        "score": frame.get("score", frame.get("source", ".")),
    }


def annotate_protein(aa: str, exons: list[dict] | None = None, v_hit: dict | None = None, i_hit: dict | None = None) -> dict:
    """Place the leader, V, I, transmembrane helix, and cytoplasmic tail."""
    seq = (aa or "").replace("*", "")
    v, i = choose_ig(seq, v_hit, i_hit)
    tms = receptor_tms(seq, v, i) if seq else []
    ig_end = (i or v)["end"] if (v or i) else 0
    dom_start = min(d["start"] for d in (v, i) if d) if (v or i) else 0
    ranged: list[dict] = []
    if exons:
        if exons[0].get("aa_start") not in (None, "", "."):
            ranged = [dict(ex) for ex in exons]
            for n, item in enumerate(ranged, start=1):
                item["exon"] = item.get("exon") or n
        else:
            for n, (ex, a1, a2) in enumerate(exon_aa_ranges(exons), start=1):
                item = dict(ex)
                item["exon"] = item.get("exon") or n
                item["aa_start"] = a1
                item["aa_end"] = a2
                ranged.append(item)
    sp, leader_flags, leader_end = leader_call(ranged, len(seq), dom_start) if ranged else (False, [], 0)
    if not sp and v and v["start"] > 12 and not ranged:
        leader_end = v["start"] - 1
        sp = True
    flags = list(leader_flags)
    if v and not v.get("g16"):
        flags.append("V missing G16")
    if i and not i.get("g16"):
        flags.append("I missing G16")
    if i and i.get("source") == "cys_frame" and (not i_hit or i.get("canonical")):
        flags.append("I placed by cysteine frame")
    if not seq.startswith("M"):
        flags.append("does not start with methionine")
    if v and i and exon_at(ranged, len(seq), v["c23"]) and exon_at(ranged, len(seq), v["c23"]) == exon_at(ranged, len(seq), i["c23"]):
        flags.append("V and I share one exon")
    return {
        "v_dom": v,
        "i_dom": i,
        "tms": tms,
        "sp": sp,
        "v": bool(v),
        "i": bool(i),
        "tm": len(tms) == 1,
        "ct": bool(tms) and tms[0][1] < len(seq),
        "leader_end": leader_end,
        "ig_end": ig_end,
        "type": signaling_type(seq, tms, ig_end) if seq else "",
        "flags": flags,
        "flag": "; ".join(flags),
        "exons": ranged,
        "seq": seq,
    }


def trim_to_leader_met(seq: str, exons: list[dict] | None = None, v_from: int = 0) -> tuple[int, str]:
    """Keep a Met already in the spliced leader exon. Never walk the genome.

    miniprot sometimes translates a few extra residues *inside the first CDS
    exon* before the real start (NSKHKDHFEGRLIMRYY… vs MRYYF…). We may clip
    those amino acids if that Met is still in exon 1, still before V, and is
    followed by a hydrophobic signal.

    We do **not** scan upstream genomic sequence for an ATG. That walk leaves
    the spliced CDS and can pull intron into the protein.
    """
    seq = (seq or "").replace("*", "")
    if not seq or seq[0] == "M":
        return 0, seq
    # Only search inside the first CDS exon, and only before the V domain.
    max_pos = min(len(seq), 30)
    if exons:
        first_nt = int(exons[0].get("seq_len") or (int(exons[0]["end"]) - int(exons[0]["start"]) + 1))
        max_pos = min(max_pos, max(1, first_nt // 3))
    if v_from and int(v_from) > 1:
        max_pos = min(max_pos, int(v_from) - 1)
    best_i = None
    best_n = -1
    for i, ch in enumerate(seq[:max_pos]):
        if ch != "M":
            continue
        n = sum(c in HYDRO for c in seq[i + 1 : i + 16])
        if n >= 6 and (n > best_n or (n == best_n and i > (best_i or -1))):
            best_n, best_i = n, i
    if best_i is not None:
        return best_i, seq[best_i:]
    return 0, seq


def trim_exons_nt(exons: list[dict], n_aa: int, strand: str) -> list[dict]:
    """Clip n_aa residues from the transcript 5' end (first CDS exon)."""
    cut = max(0, int(n_aa)) * 3
    if cut <= 0:
        return list(exons)
    kept = []
    for ex in exons:
        length = int(ex.get("seq_len") or (int(ex["end"]) - int(ex["start"]) + 1))
        if cut >= length:
            cut -= length
            continue
        ex = dict(ex)
        if cut:
            if (strand or "+") == "-":
                ex["end"] = int(ex["end"]) - cut
            else:
                ex["start"] = int(ex["start"]) + cut
            length -= cut
            cut = 0
        ex["seq_len"] = length
        kept.append(ex)
    return kept


def shift_hits(hits, trim: int) -> list[dict]:
    out = []
    for h in _as_hits(hits):
        a = int(h["ali_from"]) - trim
        b = int(h["ali_to"]) - trim
        if b < 1:
            continue
        out.append({**h, "ali_from": max(1, a), "ali_to": b})
    return out


def apply_leader_trim(row: dict) -> int:
    """Trim row['aa'] / exons / HMM coords to a Met already in exon 1. Returns aa removed."""
    aa = (row.get("aa") or "").replace("*", "")
    v_from = 0
    try:
        v_from = int(row.get("v_aa_from") or 0)
    except (TypeError, ValueError):
        v_from = 0
    trim, new = trim_to_leader_met(aa, exons=row.get("exons") or None, v_from=v_from)
    if not trim:
        row["aa"] = new
        return 0
    row["aa"] = new
    row["leader_trim"] = trim
    strand = row.get("strand") or "+"
    if row.get("exons"):
        row["exons"] = trim_exons_nt(row["exons"], trim, strand)
        if row["exons"]:
            row["start"] = min(int(e["start"]) for e in row["exons"])
            row["end"] = max(int(e["end"]) for e in row["exons"])
            row["span"] = int(row["end"]) - int(row["start"]) + 1
            row["n_exons"] = len(row["exons"])
    for key in ("v_aa_from", "v_aa_to", "i_aa_from", "i_aa_to"):
        val = row.get(key)
        if val in ("", None, "."):
            continue
        try:
            row[key] = max(1, int(val) - trim)
        except (TypeError, ValueError):
            pass
    return trim


def domain_track(aa: str, v_hit: dict | None, i_hit: dict | None) -> list[dict]:
    """Leader / V / I / TM / CYT. TM is after the last Ig, never used to block I."""
    seq = (aa or "").replace("*", "")
    if not seq:
        return []
    v, i = choose_ig(seq, v_hit, i_hit)
    tms = receptor_tms(seq, v, i)
    out: list[dict] = []
    first_ig = min(d["start"] for d in (v, i) if d) if (v or i) else 0
    if first_ig >= 12:
        out.append({"domain": "leader", "aa_start": 1, "aa_end": first_ig - 1, "score": "."})
    if v:
        out.append(
            {
                "domain": "V",
                "aa_start": v["start"],
                "aa_end": v["end"],
                "score": v.get("score", v.get("source", ".")),
            }
        )
    if i:
        out.append(
            {
                "domain": "I",
                "aa_start": i["start"],
                "aa_end": i["end"],
                "score": i.get("score", i.get("source", ".")),
            }
        )
    last_ig = max([d["aa_end"] for d in out if d["domain"] in {"V", "I"}] or [0])
    if tms:
        tm = tms[0]
        out.append({"domain": "TM", "aa_start": tm[0], "aa_end": tm[1], "score": "."})
        cyt0 = tm[1] + 1
    else:
        cyt0 = last_ig + 1 if last_ig else 0

    if cyt0 and cyt0 <= len(seq):
        tail = seq[cyt0 - 1 :]
        hits = list(ITIM.finditer(tail))
        if len(hits) >= 2:
            split = cyt0 + hits[0].end()
            split = min(max(split, cyt0 + 8), len(seq) - 8)
            out.append({"domain": "CYT1", "aa_start": cyt0, "aa_end": split, "score": "."})
            out.append({"domain": "CYT2", "aa_start": split + 1, "aa_end": len(seq), "score": "."})
        else:
            name = "CYT1" if hits or len(tail) >= 8 else "CYT"
            out.append({"domain": name, "aa_start": cyt0, "aa_end": len(seq), "score": "."})
    return out


def exon_aa_ranges(exons: list[dict]) -> list[tuple[dict, int, int]]:
    """Map CDS-order exons to protein coordinates (phase-aware via nt length)."""
    acc = 0
    ranges = []
    for ex in exons:
        nt = int(ex.get("seq_len") or (int(ex["end"]) - int(ex["start"]) + 1))
        aa_start = acc // 3 + 1
        acc += nt
        aa_end = acc // 3
        if aa_end >= aa_start:
            ranges.append((ex, aa_start, aa_end))
    return ranges


def domains_on_exon(ranges: list[tuple[dict, int, int]], domains: list[dict]) -> list[str]:
    """Labels for each exon. Ig (V/I) wins over TM on a fused exon."""
    labels = []
    for _ex, a1, a2 in ranges:
        occ = []
        for d in domains:
            n = overlap_len(a1, a2, d["aa_start"], d["aa_end"])
            if n >= 8:
                occ.append((DISPLAY.get(d["domain"], d["domain"]), n, d["domain"]))
        if not occ:
            labels.append(".")
            continue
        ig = [name for name, _n, raw in occ if raw in {"V", "I"}]
        if ig:
            labels.append("+".join(dict.fromkeys(ig)))
            continue
        labels.append("+".join(dict.fromkeys(name for name, _n, _raw in occ)))
    return labels


def parse_hmmalign_columns(sto: Path) -> dict[str, dict]:
    """G16 / C23 / Cys occupancy from HMMER match states (IMGT-style columns)."""
    seqs: dict[str, str] = {}
    rf = ""
    if not sto.exists():
        return {}
    with sto.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith("# STOCKHOLM") or line.startswith("//"):
                continue
            if line.startswith("#=GC RF"):
                rf += line.split()[-1]
                continue
            if line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 2:
                seqs[parts[0]] = seqs.get(parts[0], "") + parts[1]
    out = {}
    for name, seq in seqs.items():
        if not rf or len(seq) != len(rf):
            # RF may be shorter if the last block omitted it; trim to min
            n = min(len(seq), len(rf))
            seq, rf_use = seq[:n], rf[:n]
        else:
            rf_use = rf
        col = 0
        g16 = c23 = "."
        cys = []
        for res, mark in zip(seq, rf_use):
            if mark not in {"x", "X"}:
                continue
            col += 1
            aa = res.upper()
            if col == 16:
                g16 = aa if aa not in {".", "-", "_"} else "."
            if col == 23:
                c23 = aa if aa not in {".", "-", "_"} else "."
            if aa == "C":
                cys.append(str(col))
        out[name] = {
            "g16": g16,
            "c23": c23,
            "cys_aln_positions": ",".join(cys) if cys else ".",
            "hmm_len": col,
        }
    return out


def align_ig_columns(records: list[tuple[str, str]], hmm: Path, sto: Path) -> dict[str, dict]:
    """hmmalign peptides to a G16/C23-framed HMM; return per-id column occupancy."""
    if not records or not hmm.exists():
        return {}
    fa = sto.with_suffix(".faa")
    p.write_fasta(fa, records)
    sto.parent.mkdir(parents=True, exist_ok=True)
    try:
        with sto.open("w") as fh:
            subprocess.check_call(["hmmalign", "--amino", str(hmm), str(fa)], stdout=fh)
    except (subprocess.CalledProcessError, FileNotFoundError) as exc:
        print(f"hmmalign failed {hmm.name}: {exc}", flush=True)
        return {}
    return parse_hmmalign_columns(sto)


def sync_row_from_domains(row: dict, domains: list[dict]) -> None:
    """n_i / TM / scores come from the filtered track, not the raw HMM table."""
    i = next((d for d in domains if d["domain"] == "I"), None)
    v = next((d for d in domains if d["domain"] == "V"), None)
    row["n_i"] = 1 if i else 0
    row["n_v"] = 1 if v else int(row.get("n_v") or 0)
    row["tm"] = any(d["domain"] == "TM" for d in domains)
    if i:
        row["i_score"] = i.get("score", row.get("i_score", "."))
        row["i_aa_from"] = i["aa_start"]
        row["i_aa_to"] = i["aa_end"]
    else:
        row["i_score"] = 0
        row["i_aa_from"] = ""
        row["i_aa_to"] = ""
    if v:
        row["v_aa_from"] = v["aa_start"]
        row["v_aa_to"] = v["aa_end"]
        row["v_score"] = v.get("score", row.get("v_score", "."))


def write_exports(
    unique: list[dict],
    out: Path,
    v_hits: dict | None = None,
    i_hits: dict | None = None,
    v_hmm: Path | None = None,
    i_hmm: Path | None = None,
) -> None:
    """Write exons.tsv, domains.tsv, unique_loci.gff3, proteins/unique_complete.dom.faa."""
    tables = out / "tables"
    tables.mkdir(parents=True, exist_ok=True)
    v_hits = v_hits or {}
    i_hits = i_hits or {}

    exon_rows = []
    dom_rows = []
    gff_lines = ["##gff-version 3"]
    fa_recs = []

    for r in unique:
        gid = r.get("gene_id", ".")
        aa = (r.get("aa") or "").replace("*", "")
        vh = v_hits.get(gid) or None
        ih = i_hits.get(gid) or None
        if vh is None and r.get("v_aa_from"):
            vh = {"ali_from": r["v_aa_from"], "ali_to": r["v_aa_to"], "score": r.get("v_score", ".")}
        if ih is None and r.get("i_aa_from"):
            ih = {"ali_from": r["i_aa_from"], "ali_to": r["i_aa_to"], "score": r.get("i_score", ".")}
        vh = pick_v_hit(vh, len(aa)) if aa else None
        ih = (
            pick_i_hit(
                ih,
                int(vh["ali_from"]) if vh else 0,
                int(vh["ali_to"]) if vh else 0,
                None,
                len(aa),
                tm=None,
            )
            if aa
            else None
        )

        domains = domain_track(aa, vh, ih) if aa else []
        placed = annotate_protein(aa, r.get("exons") or [], vh, ih) if aa else {}
        sync_row_from_domains(r, domains)
        r["_domains"] = domains
        r["_placed"] = placed
        name = r.get("known") or r.get("known_overlap") or "."
        status = r.get("status", ".")
        attrs = f"ID={gid};Name={name};status={status}"
        if r.get("query"):
            attrs += f";query={r['query']}"
        gff_lines.append(
            "\t".join(
                [
                    r["contig"],
                    "nitr_catalog",
                    "gene",
                    str(r["start"]),
                    str(r["end"]),
                    ".",
                    r.get("strand") or ".",
                    ".",
                    attrs,
                ]
            )
        )

        exons = list(r.get("exons") or [])
        ranges = exon_aa_ranges(exons) if exons else []
        placed_exons = placed.get("exons") or []
        exo_labels = []
        if placed_exons and (placed.get("v_dom") or placed.get("i_dom")):
            for ex in placed_exons:
                lab, _pep = exon_display(
                    aa,
                    ex,
                    placed_exons,
                    placed.get("v_dom"),
                    placed.get("i_dom"),
                    placed.get("tms") or [],
                    placed.get("leader_end") or 0,
                )
                exo_labels.append(lab or ".")
        else:
            exo_labels = domains_on_exon(ranges, domains) if ranges else []
        if len(exo_labels) != len(ranges):
            exo_labels = domains_on_exon(ranges, domains) if ranges else []
        for i, ((ex, a1, a2), lab) in enumerate(zip(ranges, exo_labels), 1):
            phase = str(ex.get("phase", "0"))
            gff_lines.append(
                "\t".join(
                    [
                        r["contig"],
                        "nitr_catalog",
                        "CDS",
                        str(ex["start"]),
                        str(ex["end"]),
                        ".",
                        r.get("strand") or ".",
                        phase if phase in {"0", "1", "2"} else "0",
                        f"Parent={gid};exon={i};aa={a1}-{a2};domains={lab}",
                    ]
                )
            )
            exon_rows.append(
                {
                    "gene_id": gid,
                    "exon": i,
                    "contig": r["contig"],
                    "start": ex["start"],
                    "end": ex["end"],
                    "strand": r.get("strand", "."),
                    "phase": phase,
                    "aa_start": a1,
                    "aa_end": a2,
                    "domains": lab,
                    "donor": ex.get("donor", "."),
                    "acceptor": ex.get("acceptor", "."),
                    "status": status,
                    "name": name,
                }
            )

        for d in domains:
            dom_rows.append(
                {
                    "gene_id": gid,
                    "name": name,
                    "status": status,
                    "domain": d["domain"],
                    "aa_start": d["aa_start"],
                    "aa_end": d["aa_end"],
                    "score": d.get("score", "."),
                    "aa_len": len(aa),
                    "g16": ".",
                    "c23": ".",
                    "cys_aln_positions": ".",
                }
            )
        if aa:
            tag = ",".join(f"{d['domain']}:{d['aa_start']}-{d['aa_end']}" for d in domains) or "none"
            fa_recs.append((f"{gid} name={name} status={status} domains={tag}", aa))

    v_peps, i_peps = [], []
    for r in unique:
        aa = (r.get("aa") or "").replace("*", "")
        for d in r.get("_domains") or []:
            if d["domain"] not in {"V", "I"}:
                continue
            pep = aa[d["aa_start"] - 1 : d["aa_end"]]
            if len(pep) < 40:
                continue
            rec = (r["gene_id"], pep)
            if d["domain"] == "V":
                v_peps.append(rec)
            else:
                i_peps.append(rec)
    v_col = align_ig_columns(v_peps, v_hmm, tables / "V_hmmalign.sto") if v_hmm else {}
    i_col = align_ig_columns(i_peps, i_hmm, tables / "I_hmmalign.sto") if i_hmm else {}
    for row in dom_rows:
        cols = v_col if row["domain"] == "V" else i_col if row["domain"] == "I" else {}
        info = cols.get(row["gene_id"]) or {}
        if info:
            row["g16"] = info.get("g16", ".")
            row["c23"] = info.get("c23", ".")
            row["cys_aln_positions"] = info.get("cys_aln_positions", ".")
            row["v_g16c23_hmm" if row["domain"] == "V" else "i_g16c23_hmm"] = (
                info.get("g16") == "G" and info.get("c23") == "C"
            )
            if row["domain"] == "V":
                rec = next((u for u in unique if u.get("gene_id") == row["gene_id"]), None)
                if rec:
                    rec["v_g16c23_hmm"] = info.get("g16") == "G" and info.get("c23") == "C"
            if row["domain"] == "I":
                rec = next((u for u in unique if u.get("gene_id") == row["gene_id"]), None)
                if rec:
                    rec["i_g16c23_hmm"] = info.get("g16") == "G" and info.get("c23") == "C"

    def dump(path: Path, rows: list[dict], cols: list[str]) -> None:
        with path.open("w") as fh:
            fh.write("\t".join(cols) + "\n")
            for row in rows:
                fh.write("\t".join(str(row.get(c, ".")) for c in cols) + "\n")

    dump(
        tables / "exons.tsv",
        exon_rows,
        [
            "gene_id",
            "name",
            "status",
            "exon",
            "contig",
            "start",
            "end",
            "strand",
            "phase",
            "aa_start",
            "aa_end",
            "domains",
            "donor",
            "acceptor",
        ],
    )
    dump(
        tables / "domains.tsv",
        dom_rows,
        [
            "gene_id",
            "name",
            "status",
            "domain",
            "aa_start",
            "aa_end",
            "score",
            "aa_len",
            "g16",
            "c23",
            "cys_aln_positions",
        ],
    )
    (tables / "unique_loci.gff3").write_text("\n".join(gff_lines) + "\n", encoding="utf-8")
    p.write_fasta(out / "proteins" / "unique_complete.dom.faa", fa_recs)
    print(
        f"Wrote exons={len(exon_rows)} domains={len(dom_rows)} gff genes={sum(1 for r in unique)}",
        flush=True,
    )
