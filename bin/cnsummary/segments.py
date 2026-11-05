"""Segments: reading them and preparing them for the calls.

A segment is a dict. After prepare() every segment has these keys:

  chrom, start, end      position (1-based, inclusive, no "chr")
  major, minor           the tool's integer state
  a_raw, b_raw           unrounded copies per tumour cell, a_raw >= b_raw (the integer state without segments_raw)
  n_het                  germline het SNPs in the segment, or None
  a, b                   a_raw, b_raw with the allele scale divided out (see estimate_scale)
  c_major, c_minor       the state after the balance correction (see correct_balanced); every call uses these
  lowhet                 het-SNP density far below the median: the normal is homozygous too (germline, not LOH)
  masked                 >= 50% inside a --mask region
"""
from collections import defaultdict

from . import settings
from .genome import AUTOSOMES, ORDER, length, norm_chrom
from .tables import read_tsv


def load(path, have_raw):
    """Segments from <prefix>.segments_raw.tsv (have_raw) or <prefix>.segments.tsv, sorted by position."""
    segs = []
    for r in read_tsv(path):
        ma, mi = int(float(r["nMajor"])), int(float(r["nMinor"]))
        a, b = (float(r["nAraw"]), float(r["nBraw"])) if have_raw else (float(ma), float(mi))
        nh = r.get("n_het")
        # positions via float: R writes round numbers like 8000000 as 8e+06
        segs.append({"chrom": norm_chrom(r["chr"]), "start": int(float(r["startpos"])), "end": int(float(r["endpos"])),
                     "major": ma, "minor": mi, "a_raw": max(a, b), "b_raw": min(a, b),
                     "n_het": float(nh) if nh not in (None, "", "NA") else None})
    return sorted(segs, key=lambda d: (ORDER.get(d["chrom"], 99), d["start"]))


def estimate_scale(segs, min_len=5000000):
    """When purity is under-fitted, every allele's distance from 1 copy is inflated by the same factor. It is the
    median (major - 1) over 2+1 segments >= min_len: 1.0 for an unbiased fit, capped at SCALE_MAX."""
    xs = sorted(d["a_raw"] - 1 for d in segs if (d["major"], d["minor"]) == (2, 1) and length(d) >= min_len)
    if len(xs) < settings.SCALE_MIN_SEGMENTS:
        return 1.0
    med = xs[len(xs) // 2] if len(xs) % 2 else (xs[len(xs) // 2 - 1] + xs[len(xs) // 2]) / 2
    return min(max(med, 1.0), settings.SCALE_MAX)


def correct_balanced(segs, tol):
    """An allelically balanced segment (|a - b| < tol) gets an even state round(mean)+round(mean): an odd state
    there is a rounding artefact of a low purity fit (e.g. tetrasomy as 3+2)."""
    for d in segs:
        d["c_major"], d["c_minor"] = d["major"], d["minor"]
        if d["major"] != d["minor"] and abs(d["a"] - d["b"]) < tol and d["b"] >= 0.5:
            k = max(1, round((d["a"] + d["b"]) / 2))
            d["c_major"] = d["c_minor"] = k


def flag_artefacts(segs, mask, density_min):
    """Sets lowhet (het-SNP density < density_min x the median over autosomal segments >= 1 Mb) and masked
    (>= 50% inside a mask region). Returns whether the density filter could be applied (n_het present)."""
    have = all(d["n_het"] is not None for d in segs) and segs
    dens = sorted(d["n_het"] / (length(d) / 1e6) for d in segs
                  if have and d["chrom"] in AUTOSOMES and length(d) >= 1000000)
    med = dens[len(dens) // 2] if dens else 0.0
    for d in segs:
        n = length(d)
        d["lowhet"] = bool(have and med > 0 and d["n_het"] / (n / 1e6) < density_min * med)
        d["masked"] = any(d["chrom"] == c and min(d["end"], e) - max(d["start"], s) >= 0.5 * n for c, s, e in mask)
    return bool(have)


def prepare(segs, have_raw, args, mask):
    """Every preparation step, in order: allele scale, balance correction, artefact flags.
    -> (scale, have_nhet)."""
    scale = estimate_scale(segs) if have_raw else 1.0
    for d in segs:
        d["a"] = 1 + (d["a_raw"] - 1) / scale
        d["b"] = 1 + (d["b_raw"] - 1) / scale
    correct_balanced(segs, args.balanced_tol if have_raw else -1)
    have_nhet = flag_artefacts(segs, mask, args.het_density_min)
    return scale, have_nhet


def merge(items, max_gap, key):
    """Merges consecutive items on the same chromosome with the same key() and a gap <= max_gap. Numeric values
    (a, b, a_raw, b_raw, cn, frac, bal) are length-weighted."""
    out = []
    for d in items:
        prev = out[-1] if out else None
        if prev and prev["chrom"] == d["chrom"] and key(prev) == key(d) and d["start"] - prev["end"] <= max_gap:
            n0, n1 = length(prev), length(d)
            for k in ("a", "b", "a_raw", "b_raw", "cn", "frac", "bal"):
                if k in d:
                    prev[k] = (prev[k] * n0 + d[k] * n1) / (n0 + n1)
            prev["end"] = d["end"]
        else:
            out.append(dict(d))
    return out


def per_chrom(segs, major, minor):
    """{chrom: {cn, loh, state, state_frac}}: length-weighted mean copies, LOH fraction and the modal state, using
    the given state keys ("major"/"minor": the tool's; "c_major"/"c_minor": corrected)."""
    L, W, Lloh = defaultdict(int), defaultdict(float), defaultdict(int)
    states = defaultdict(lambda: defaultdict(int))
    for d in segs:
        c, n = d["chrom"], length(d)
        ma, mi = d[major], d[minor]
        L[c] += n
        W[c] += (ma + mi) * n
        if mi == 0:
            Lloh[c] += n
        states[c]["%d+%d" % (ma, mi)] += n
    out = {}
    for c in L:
        state, n = max(states[c].items(), key=lambda kv: kv[1])
        out[c] = {"cn": W[c] / L[c], "loh": Lloh[c] / L[c], "state": state, "state_frac": n / L[c]}
    return out
