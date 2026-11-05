"""The standardised layers of cn_summary.json: one row per chromosome arm, and the reportable events."""
from collections import defaultdict

from . import settings
from .genome import ACROCENTRIC, AUTOSOMES, CENTROMERE, arm_part, expected_cn, length, overlap


def arm_table(segs, sex, arms):
    """One row per chromosome arm (acrocentrics: q only), the same rows for every sample and tool.
    class, relative to the expected copies (2; 1 for a male X): normal / gain / loss / cnLOH / LOH+gain /
    amp (>= expected + 3) / homdel / mixed (no total over >= half the arm). LOH on low-het or masked segments is
    not counted. major_raw/minor_raw are length-weighted over the segments with the arm's modal total, so a
    partial-arm event keeps its own cell fraction."""
    rows = []
    for c in AUTOSOMES + ["X"]:
        for arm in ("p", "q"):
            if arm == "p" and c in ACROCENTRIC:
                continue
            parts = []
            for d in segs:
                if d["chrom"] != c:
                    continue
                s, e = arm_part(c, arm, d["start"], d["end"])
                if e > s:
                    parts.append((d, e - s, s, e))
            if not parts:
                continue
            L = sum(n for _, n, _, _ in parts)
            exp = expected_cn(c, sex)
            male_x = c == "X" and sex == "XY"
            loh = 0 if male_x else sum(n for d, n, _, _ in parts
                                       if d["c_minor"] == 0 and not d.get("lowhet") and not d.get("masked"))
            total, frac, state = arms[(c, arm)]
            minor = int(state.split("+")[1])
            lohish = minor == 0 and not male_x and loh / L >= settings.ARM_MIN
            if frac < settings.ARM_MIN:
                cls = "mixed"
            elif total == 0:
                cls = "homdel"
            elif total < exp:
                cls = "loss"
            elif total >= exp + 3:
                cls = "amp"
            elif total > exp:
                cls = "LOH+gain" if lohish else "gain"
            elif lohish:
                cls = "cnLOH"
            else:
                cls = "normal"
            modal = [p for p in parts if p[0]["c_major"] + p[0]["c_minor"] == total] or parts
            Lm = sum(n for _, n, _, _ in modal)
            rows.append(dict(chrom=c, arm=arm, start=min(s for _, _, s, _ in parts), end=max(e for _, _, _, e in parts),
                             expected=exp, state=state, total=total, arm_fraction=round(frac, 3), cls=cls,
                             loh_fraction=round(loh / L, 3),
                             modal_start=min(s for _, _, s, _ in modal), modal_end=max(e for _, _, _, e in modal),
                             major_raw=round(sum(d["a"] * n for d, n, _, _ in modal) / Lm, 3),
                             minor_raw=round(sum(d["b"] * n for d, n, _, _ in modal) / Lm, 3),
                             n_het=int(sum((d["n_het"] or 0) * n / max(length(d), 1) for d, n, _, _ in parts))))
    return rows


def arm_fraction_of_cells(r):
    """Fraction of tumour cells with an arm's change, from its unrounded copies vs its integer state."""
    exp, tot = r["expected"], r["major_raw"] + r["minor_raw"]
    if r["cls"] in ("gain", "LOH+gain", "amp") and r["total"] > exp:
        f = (tot - exp) / (r["total"] - exp)
    elif r["cls"] in ("loss", "homdel") and r["total"] < exp:
        f = (exp - tot) / (exp - r["total"])
    elif r["cls"] == "cnLOH":
        f = 1 - r["minor_raw"]
    else:
        return None
    return min(max(f, 0.0), 1.0)


def arms_of(chrom, start, end, min_overlap=1000000, min_share=0.1):
    """Arm names ("7p,7q") an interval overlaps by >= min_overlap and >= min_share of its length, so an event that
    starts just before the centromere is not labelled with both arms."""
    if chrom not in CENTROMERE:
        return ""
    out = []
    for arm in ("p", "q"):
        s, e = arm_part(chrom, arm, start, end)
        if not (arm == "p" and chrom in ACROCENTRIC) and e - s >= max(min_overlap, min_share * (end - start + 1)):
            out.append(chrom + arm)
    return ",".join(out)


SUBCLONAL_TYPE = {"subclonal cnLOH": "cnLOH", "subclonal LOH": "LOH+loss", "subclonal gain": "gain"}


def build_events(arm_rows, loh, sub, clonal_min, share, tcf, subclonal_max=0.7, seg_gains=()):
    """Reportable events. Chromosome scope when all arms (acrocentric: q) share one changed class and total, else
    one arm event per changed arm; then the LOH segments, segmental gains and subclonal events (segment scope)
    that an arm/chromosome event of the same type does not already cover (>= COVERED of their length).
    clonality: clonal (>= clonal_min of tumour cells), subclonal (< subclonal_max), possibly_subclonal in between,
    since arm fractions scatter by ~0.1 and one cut-off would flip labels on noise.
    share(f), tcf(f): the share of all cells and the tumour-cell fraction to report for a fraction f."""
    def clonality(f):
        if f is None or f >= clonal_min:
            return "clonal"
        return "possibly_subclonal" if f >= subclonal_max else "subclonal"
    events = []
    by_chrom = defaultdict(list)
    for r in arm_rows:
        by_chrom[r["chrom"]].append(r)
    for c in AUTOSOMES + ["X"]:
        rs = by_chrom.get(c, [])
        changed = [r for r in rs if r["cls"] not in ("normal", "mixed")]
        if not changed:
            continue
        whole = len(changed) == len(rs) and len({r["cls"] for r in rs}) == 1 and len({r["total"] for r in rs}) == 1
        for g in ([rs] if whole else [[r] for r in changed]):
            fs = [x for x in (arm_fraction_of_cells(r) for r in g) if x is not None]
            f = sum(fs) / len(fs) if fs else None
            events.append(dict(scope="chromosome" if whole else "arm", type=g[0]["cls"], chrom=c,
                               start=min(r["modal_start"] for r in g), end=max(r["modal_end"] for r in g),
                               arms=",".join(c + r["arm"] for r in g), state=g[0]["state"], copies=g[0]["total"],
                               major_raw=g[0]["major_raw"] if len(g) == 1 else round(sum(r["major_raw"] for r in g) / len(g), 3),
                               minor_raw=g[0]["minor_raw"] if len(g) == 1 else round(sum(r["minor_raw"] for r in g) / len(g), 3),
                               clonality=clonality(f), tumour_cell_fraction=tcf(f), cell_share=share(f)))
    big = list(events)

    def covered(typ, d):
        return any(e["type"] == typ and overlap(e, d) >= settings.COVERED * length(d) for e in big)

    def segment(d, typ, state, copies, a, b, f):
        return dict(scope="segment", type=typ, chrom=d["chrom"], start=d["start"], end=d["end"],
                    arms=arms_of(d["chrom"], d["start"], d["end"]), state=state, copies=copies,
                    major_raw=None if a is None else round(a, 3), minor_raw=None if b is None else round(b, 3),
                    clonality=clonality(f), tumour_cell_fraction=tcf(f), cell_share=share(f))
    for d in loh:
        if not covered(d["type"], d):
            events.append(segment(d, d["type"], "", int(round(d["cn"])), None, None, d["frac"]))
    for d in seg_gains:
        if not covered("gain", d):
            events.append(segment(d, "gain", d["state"], int(round(d["total"])), d["a"], d["b"], d["frac"]))
    for d in sub:
        typ = SUBCLONAL_TYPE.get(d["kind"], d["kind"])
        if not covered(typ, d):
            events.append(segment(d, typ, d["state"], None, d["a"], d["b"], min(d["cn"], 1.0)))
    for e in events:
        e["length_mb"] = round(length(e) / 1e6, 2)
        for k in ("tumour_cell_fraction", "cell_share"):
            if e[k] is not None:
                e[k] = round(e[k], 3)
    return events
