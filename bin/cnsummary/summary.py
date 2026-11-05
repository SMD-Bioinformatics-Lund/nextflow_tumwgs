"""summarise(): from the fit and the segments to every call, in one dict that report.py writes out."""
import os

from . import calls, count, qc, segments
from .genome import AUTOSOMES, parse_regions
from .tables import num, read_tsv


def summarise(args):
    """args: the parsed command line (see cn_summary.py). -> dict s with every result; s["segs"] empty = no fit."""
    s = {}
    fit = s["fit"] = read_tsv(args.purity_ploidy)[0]
    s["gate"] = read_tsv(args.gate)[0] if args.gate and os.path.exists(args.gate) else {}
    have_raw = s["have_raw"] = (bool(args.segments_raw) and os.path.exists(args.segments_raw)
                                and os.path.getsize(args.segments_raw) > 0)
    segs = s["segs"] = segments.load(args.segments_raw if have_raw else args.segments, have_raw)
    s["mask"] = parse_regions(args.mask)
    s["scale"], s["have_nhet"] = segments.prepare(segs, have_raw, args, s["mask"])
    s["germline_roh"] = calls.germline_homozygosity(segs, args.sex, args.min_loh_len, args.max_gap)

    # purity: reported only when measured; otherwise fractions are shares of all cells (fit purity x fraction
    # does not depend on the arbitrary fit purity)
    pur = s["purity"] = num(fit.get("purity"))
    sub = calls.subclonal_segments(segs, args.subclonal_min_dev, args.subclonal_min_len, args.max_gap,
                                   args.subclonal_min_dev_loh, args.balanced_gain_min_len) if have_raw else []
    loh = s["loh"] = calls.loh_segments(segs, args.sex, args.min_loh_len, args.max_gap)
    s["purity_status"], s["purity_min"] = qc.purity_state(fit.get("purity_status") or "ok", pur, segs, loh, sub,
                                                          args.sex, args.subclonal_min_len)
    pur_rep = s["purity_reported"] = pur if s["purity_status"] in ("ok", "manual") else None
    all_cells = pur_rep is None
    s["share"] = (lambda f: f * pur) if all_cells and pur else (lambda f: f)
    s["cells"] = "all cells" if all_cells else "tumour cells"
    s["cell_share"] = lambda f: None if f is None or pur is None else f * pur
    s["tcf"] = lambda f: None if f is None or pur_rep is None else f
    s["per_chrom_tool"] = segments.per_chrom(segs, "major", "minor")
    pc = s["per_chrom"] = segments.per_chrom(segs, "c_major", "c_minor")

    s["seg_gains"] = []
    if segs:
        arms = s["arms"] = calls.per_arm(segs, "c_major", "c_minor")
        copies = s["copies"] = count.chromosome_copies(arms, args.sex, args.whole_frac)
        s["seg_gains"] = calls.segmental_gains(segs, args.sex, calls.gained_arms(arms, copies, args.sex),
                                               args.subclonal_min_len, args.max_gap)
        sub = calls.drop_inside(sub, s["seg_gains"])
    s["sub"] = sub

    # QC of the pair
    s["floor"] = num(fit.get("purity_floor"))
    s["xr"] = qc.x_het_ratio(segs) if s["have_nhet"] else None
    s["sex_data"], s["sex_status"] = qc.sex_check(fit, args.sex, s["xr"])
    s["tin"] = num(fit.get("tin_fraction"))
    ys = s["y"] = qc.y_status(args.ydepth, segs, pur_rep) if args.sex == "XY" else None
    s["y_present"] = not (ys and ys["lost"] >= 0.5)
    s["y_border"] = bool(ys and abs(ys["lost"] - 0.5) < args.borderline_margin)
    s["corrections"] = [d for d in segs if (d["major"], d["minor"]) != (d["c_major"], d["c_minor"])]

    # chromosome count and verdict
    s.update(count=None, count_tool=None, borderline={}, range=None, purity_alt=[], key_tri={}, triple=False,
             gained=[], lost=[], verdict_code="not_assessable", verdict="no segments")
    if segs:
        def whole(c):
            return c in pc and pc[c]["state_frac"] >= args.whole_frac

        def balanced(c):
            ma, mi = pc[c]["state"].split("+")
            return ma == mi
        copies = s["copies"]
        s["whole_loh"] = [c for c in AUTOSOMES if whole(c) and pc[c]["state"].endswith("+0")]
        s["gained"] = [c for c in AUTOSOMES if copies[c] >= 3]
        s["lost"] = [c for c in AUTOSOMES if copies[c] <= 1]
        s["gain_imbal"] = [c for c in s["gained"] if whole(c) and not balanced(c)]
        s["gain_bal"] = [c for c in s["gained"] if whole(c) and balanced(c)]
        n = s["count"] = count.chromosome_count(copies, args.sex, s["y_present"])
        copies_tool = count.chromosome_copies(calls.per_arm(segs, "major", "minor"), args.sex, args.whole_frac)
        s["count_tool"] = count.chromosome_count(copies_tool, args.sex, s["y_present"])
        s["arm_level"] = count.arm_changes(s["arms"], copies)
        s["borderline"] = count.borderline_chromosomes(s["arms"], copies, args.sex, args.whole_frac,
                                                       args.borderline_margin)
        s["purity_alt"] = (count.purity_dependent_gains(copies, sub, segs, args.sex, args.whole_frac)
                           if s["purity_status"] in ("lower_bound", "not_identifiable") else [])
        s["range"] = count.count_range(n, s["count_tool"], s["borderline"], s["purity_alt"], copies, args.sex,
                                       s["y_present"], s["y_border"])
        s["key_tri"], s["triple"] = count.heh_trisomies(copies)
        s["verdict_code"], s["verdict"] = count.verdict(n, s["whole_loh"], args.heh_verdict, args.masked_hypo_min_loh)
    return s
