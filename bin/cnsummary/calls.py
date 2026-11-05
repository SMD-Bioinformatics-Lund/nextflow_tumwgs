"""Segment-level calls: LOH, possible germline homozygosity, subclonal events and segmental gains, plus the
per-arm modal states they are checked against."""
from collections import defaultdict

from . import settings
from .genome import ACROCENTRIC, CENTROMERE, arm_part, expected_cn, length, overlap
from .segments import merge


def loh_segments(segs, sex, min_len, max_gap):
    """LOH: integer minor = 0, not low-het and not masked; a single X in a male is not LOH.
    -> [dict(chrom, start, end, cn, frac, type)], frac = share of tumour cells with the LOH (1 - remaining minor).
    Pieces shorter than max_gap are dropped first, so a short blip does not split an event in two; then only
    neighbours with the same copy number merge (a 1+0 loss next to a 2+0 cnLOH stays two events)."""
    loh = [dict(chrom=d["chrom"], start=d["start"], end=d["end"], cn=float(d["c_major"] + d["c_minor"]),
                frac=min(max(1 - d["b"], 0.0), 1.0))
           for d in segs if d["c_minor"] == 0 and not d.get("lowhet") and not d.get("masked") and length(d) >= max_gap]
    loh = merge(loh, max_gap, key=lambda d: d["cn"])
    for d in loh:
        d["type"] = "LOH+loss" if d["cn"] < 1.5 else ("cnLOH" if d["cn"] <= 2.5 else "LOH+gain")
    return [d for d in loh if length(d) >= min_len and not (sex == "XY" and d["chrom"] == "X")]


def germline_homozygosity(segs, sex, min_len, max_gap):
    """LOH-like segments with too few het SNPs: the normal is homozygous there too, so the tumour lost nothing."""
    return merge([dict(chrom=d["chrom"], start=d["start"], end=d["end"]) for d in segs
                  if d["lowhet"] and d["c_minor"] == 0 and length(d) >= min_len
                  and not (sex == "XY" and d["chrom"] == "X")], max_gap, key=lambda d: True)


def subclonal_segments(segs, min_dev, min_len, max_gap, min_dev_loh=None, balanced_min_len=20000000):
    """Segments whose alleles sit between integers, on a background with one minor copy (A+1).

    With the change in a fraction f of tumour cells:
      subclonal cnLOH : minor down by f, major up by f        (total unchanged)
      subclonal LOH   : minor down by f, major unchanged
      subclonal gain  : major up by f, minor unchanged; or both alleles up by f/2 (balanced gain k+k)
    Segments with integer minor = 0 are LOH and reported by loh_segments.
    min_dev applies to gains (read-depth driven, noisier), min_dev_loh to a lost minor allele (allele-ratio driven).
    A balanced gain has no allelic shift, so only read depth supports it: it needs >= balanced_min_len of balanced
    gain on its chromosome, which removes 5-10 Mb depth-noise pieces and keeps whole-chromosome gains.
    -> [dict(chrom, start, end, kind, state, a, b, a_raw, b_raw, cn, balanced)], cn = fraction of tumour cells."""
    loh_dev = min_dev if min_dev_loh is None else min_dev_loh
    flagged = []
    for d in segs:
        a, b = d["a"], d["b"]
        if d["c_minor"] == 0 or not 0.2 < b < 1.2 or d.get("masked"):
            continue
        db = max(1 - b, 0.0)                    # minor allele lost in a fraction db of cells
        base = max(1, int(a - db + 0.25))       # background major
        da = a - base                           # major allele gained in a fraction da of cells
        state = "%d+1" % base
        k = max(1, round(b))                    # balanced background k+k
        extra = a + b - 2 * k                   # total copies above it
        if db >= loh_dev and d.get("lowhet"):
            continue                             # allelic imbalance without het SNPs: not tumour LOH
        balanced = False
        if db >= loh_dev and abs(da - db) < max(min_dev, loh_dev):
            kind, frac = "subclonal cnLOH", (da + db) / 2
        elif db >= loh_dev and da < min_dev / 2:
            kind, frac = "subclonal LOH", db
        elif abs(a - b) < min_dev and min_dev <= extra <= 1 - min_dev:
            # a gain split evenly over both alleles (a segmenter can set a small allele shift to exactly 0.5);
            # tested before the one-allele gain, as 1.16+1.16 is a gain of 0.32, not 0.16 on the major allele
            kind, frac, state, balanced = "subclonal gain", extra, "%d+%d" % (k, k), True
        elif min_dev <= da <= 1 - min_dev and db < loh_dev:
            # a minor allele slightly below 1 (< loh_dev) is too little loss to call and must not split a gain
            kind, frac = "subclonal gain", da
        else:
            continue
        flagged.append(dict(chrom=d["chrom"], start=d["start"], end=d["end"], kind=kind, state=state,
                            a=a, b=b, a_raw=d["a_raw"], b_raw=d["b_raw"], cn=frac, bal=float(balanced)))
    flagged = merge(flagged, max_gap, key=lambda d: (d["kind"], d["state"]))
    flagged = [d for d in flagged if length(d) >= min_len]
    for d in flagged:   # merged: balanced when most of its length is
        d["balanced"] = d.pop("bal") >= 0.5
    bal = defaultdict(int)
    for d in flagged:
        if d["balanced"]:
            bal[d["chrom"]] += length(d)
    return [d for d in flagged if not (d["balanced"] and bal[d["chrom"]] < balanced_min_len)]


def per_arm(segs, major, minor):
    """{(chrom, "p"|"q"): (modal total, fraction of the arm it covers, its most common state "a+b")}."""
    L = defaultdict(lambda: defaultdict(int))
    st = defaultdict(lambda: defaultdict(int))
    for d in segs:
        c = d["chrom"]
        if c not in CENTROMERE:
            continue
        for arm in ("p", "q"):
            s, e = arm_part(c, arm, d["start"], d["end"])
            if e > s:
                L[(c, arm)][d[major] + d[minor]] += e - s
                st[(c, arm)]["%d+%d" % (d[major], d[minor])] += e - s
    out = {}
    for k, tot in L.items():
        n = sum(tot.values())
        cn, m = max(tot.items(), key=lambda kv: kv[1])
        state = max(((s, v) for s, v in st[k].items() if sum(map(int, s.split("+"))) == cn), key=lambda kv: kv[1])[0]
        out[k] = (cn, m / n, state)
    return out


def gained_arms(arms, copies, sex):
    """Arms already reported as gained: modal total above the expected copies, over >= half the arm or on a
    chromosome counted as gained."""
    return {(c, arm) for (c, arm), v in arms.items() if v[0] > expected_cn(c, sex)
            and (v[1] >= settings.ARM_MIN or copies[c] > expected_cn(c, sex))}


def segmental_gains(segs, sex, gained, min_len, max_gap):
    """Gains with an integer state above the expected copies that are neither an arm nor a chromosome gain
    (a clonal gain of less than half an arm). Pieces the subclonal classifier also flags are included, so one
    gain is one event. gained: arms already reported as gained (gained_arms); a gain with < min_len outside them
    (and outside the centromere) is not repeated.
    -> [dict(chrom, start, end, state, total, a, b, frac)], frac = tumour cells with the gain, from the raw total."""
    out = []
    for d in segs:
        if d["chrom"] not in CENTROMERE or d.get("masked"):
            continue
        exp, tot = expected_cn(d["chrom"], sex), d["c_major"] + d["c_minor"]
        if tot <= exp:
            continue
        n = length(d)
        st = "%d+%d" % (d["c_major"], d["c_minor"])
        prev = out[-1] if out else None
        if prev and prev["chrom"] == d["chrom"] and d["start"] - prev["end"] <= max_gap:
            n0 = length(prev)
            for k, v in (("a", d["a"]), ("b", d["b"]), ("total", tot)):
                prev[k] = (prev[k] * n0 + v * n) / (n0 + n)
            prev["states"][st] += n
            prev["end"] = d["end"]
        else:
            out.append(dict(chrom=d["chrom"], start=d["start"], end=d["end"], a=d["a"], b=d["b"], total=float(tot),
                            states=defaultdict(int, {st: n})))
    keep = []
    for d in out:
        outside = sum(max(0, e - s) for arm in ("p", "q") if (d["chrom"], arm) not in gained
                      and not (arm == "p" and d["chrom"] in ACROCENTRIC)
                      for s, e in [arm_part(d["chrom"], arm, d["start"], d["end"])])
        if length(d) < min_len or outside < min_len:
            continue
        exp = expected_cn(d["chrom"], sex)
        d["state"] = max(d["states"].items(), key=lambda kv: kv[1])[0]
        d["frac"] = min(max((d["a"] + d["b"] - exp) / (d["total"] - exp), 0.0), 1.0)
        del d["states"]
        keep.append(d)
    return keep


def drop_inside(sub, seg_gains):
    """Subclonal gains with >= half their length inside a segmental gain: reported there, once."""
    return [d for d in sub if not (d["kind"] == "subclonal gain" and
                                   any(overlap(g, d) >= 0.5 * length(d) for g in seg_gains))]
