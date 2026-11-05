"""Writing a summary out: the report text (heh.txt), the per-call TSVs and the standardised cn_summary.json with
its arms.tsv and events.tsv. Nothing is computed here beyond formatting."""
import json
import os

from . import qc, settings
from .count import COUNT_BASIS, count_class
from .genome import CHROMS, ORDER, length
from .layers import arm_table, build_events
from .tables import num

SCHEMA_VERSION = "1.3"


def _num4(x):
    v = num(x)
    return None if v is None else round(v, 4)


def _regions(items, label):
    """"; chrom:start-end <label(d)> ..." for a report line."""
    return "".join("; %s:%d-%d %s" % (d["chrom"], d["start"], d["end"], label(d)) for d in items)


def text_lines(s, args):
    """The report items, one "key: value" string each; "; " separates the parts of a long value (see report_text)."""
    fit, gate, pstatus, pur = s["fit"], s["gate"], s["purity_status"], s["purity"]
    share, cells = s["share"], s["cells"]
    lines = ["%s fit: purity=%s ploidy=%s goodness_of_fit=%s" % (
        args.label.upper(),
        fit.get("purity") if s["purity_reported"] is not None else "NA (%s; fit value %s)" % (pstatus, fit.get("purity")),
        fit.get("ploidy"), fit.get("goodness_of_fit"))]
    if fit.get("note"):
        lines.append("    " + fit["note"])
    if gate:
        lines.append("refit_gate: %s — %s; SNV purity %s (%s-%s, %s SNVs)" % (
            gate["decision"], gate["reason"], gate["snv_purity"], gate["snv_purity_lo"], gate["snv_purity_hi"],
            gate["n_snv"]))
    lines.append("purity_status: %s%s" % (pstatus, "" if s["purity_min"] is None
                                          else " (purity >= %.2f)" % s["purity_min"]))
    if s["floor"] is not None and pur is not None:
        lines.append("purity_floor: %.2f (lowest purity at which no confident segment needs negative copies)%s" % (
            s["floor"], " — FIT %.2f IS BELOW IT: a clonal deletion says purity is higher" % pur
            if pur < s["floor"] - settings.FLOOR_MARGIN else ""))
    if pstatus == "not_identifiable":
        lines.append("note: purity not identifiable (too little aberrant genome) — event fractions below are "
                     "shares of all cells")
    elif pstatus == "lower_bound":
        lines.append("note: gains only, no loss or LOH — purity is only bounded (>= %.2f); a clonal gain at lower "
                     "purity fits as well as a subclonal gain at purity 1; event fractions below are shares of all "
                     "cells" % s["purity_min"])
    if not s["have_nhet"]:
        lines.append("note: no n_het in the segments — germline-homozygosity filter not applied (rerun the fit)")
    lines.append("sex_check: given %s; %s%s" % (
        args.sex, {"ok": "data agree",
                   "mismatch": "DATA SAY %s — X copies, count and HeH are unreliable" % (s["sex_data"] or "XY"),
                   "unclear": "data unclear"}.get(s["sex_status"], "not assessed"),
        "" if s["xr"] is None else " (X/autosome het-SNP density %.2f)" % s["xr"]))
    if s["tin"] is not None:
        lines.append("tumour_in_normal: ~%.0f%% of normal cells%s%s" % (
            100 * s["tin"], " (se %s)" % fit.get("tin_se") if fit.get("tin_se") not in (None, "", "NA") else "",
            " — LOH and copy changes are understated; low-level events may be missed"
            if s["tin"] >= args.tin_max else ""))
    ys = s["y"]
    if args.sex == "XY":
        lines.append("loss_of_y: %s" % ("not assessed (no Y read depth)" if ys is None else
                     "Y/autosome depth %.2fx the normal; %s%s" % (
                         ys["ratio"], "Y lost in ~%.0f%% of %s" % (100 * ys["lost"], ys["cells"]) if ys["called"]
                         else "no loss (< %.0f%% of cells)" % (100 * settings.Y_LOSS_MIN),
                         "" if ys["y_copies"] is None else " (%.2f Y copies per tumour cell)" % ys["y_copies"])))
    if not s["segs"]:
        lines.append("verdict: NO FIT — no segments; HeH/LOH not assessed")
        return lines

    ploidy = num(fit.get("ploidy"))
    if ploidy is not None and ploidy > settings.HIGH_PLOIDY:
        lines.append("note: ploidy > 3.2 — may be the doubled solution of a near-diploid genome; check the sunrise plot")
    if s["scale"] > settings.SCALE_FLAG:
        lines.append("note: gained alleles read %.2fx high (2+1 major = %.2f) — purity is probably under-fitted; "
                     "raw values rescaled before subclonal/balance calls" % (s["scale"], 1 + s["scale"]))
    n, n_tool, lo_hi = s["count"], s["count_tool"], s["range"]
    lines += [
        "chromosome_count: %d (%s; copy-number based, not a karyotype)%s%s" % (
            n, count_class(n), "" if n == n_tool else "; %d from the tool's integer states" % n_tool,
            "; Y assumed present, not assessed" if args.sex == "XY" and ys is None else
            "; Y not counted (lost in most tumour cells)" if not s["y_present"] else ""),
        "count_range: %s" % ("-" if lo_hi[0] == lo_hi[1] else "%d-%d; %s" % (
            lo_hi[0], lo_hi[1],
            "; ".join(["%s counted %d, could be %d" % (c, cnt, alt) for c, (cnt, alt) in s["borderline"].items()] +
                      (["%d before the balanced correction" % n_tool] if n != n_tool else []) +
                      (["Y lost in ~half of the tumour cells"] if s["y_border"] else []) +
                      (["purity only bounded: +%s could be whole or subclonal" % ",".join(c for c, _ in s["purity_alt"])]
                       if s["purity_alt"] else []))))]
    if args.heh_verdict == "yes":
        lines.append("heh_key_trisomies: %s%s" % (
            " ".join(("+" if g else "no +") + c for c, g in s["key_tri"].items()),
            "; triple trisomy (+4, +10, +17)" if s["triple"] else ""))
    corr = s["corrections"]
    lines += [
        "balanced_corrections: %s" % ("-" if not corr else "%d segments; %s" % (
            len(corr), ", ".join("%s:%d-%d %d+%d->%d+%d" % (d["chrom"], d["start"], d["end"], d["major"], d["minor"],
                                                           d["c_major"], d["c_minor"]) for d in corr))),
        "gained_autosomes: %s" % (",".join(s["gained"]) or "-"),
        "lost_autosomes: %s" % (",".join(s["lost"]) or "-"),
        "arm_level_changes: %s (not counted)" % (", ".join(s["arm_level"]) or "-"),
        "gained_whole_imbalanced: %s" % (",".join(s["gain_imbal"]) or "-"),
        "gained_whole_balanced: %s" % (",".join(s["gain_bal"]) or "-"),
        "whole_autosome_loh: %s" % (",".join(s["whole_loh"]) or "-"),
        "segmental_gains: %d (>= %.0f Mb, part of an arm; not counted)%s" % (
            len(s["seg_gains"]), args.subclonal_min_len / 1e6,
            _regions(s["seg_gains"], lambda d: "%s in ~%.0f%% of %s" % (d["state"], 100 * share(d["frac"]), cells))),
        "loh_segments: %d (nMinor = 0, >= %.1f Mb)%s" % (
            len(s["loh"]), args.min_loh_len / 1e6,
            _regions(s["loh"], lambda d: "%s in ~%.0f%% of %s" % (d["type"], 100 * share(d["frac"]), cells))),
        "possible_germline_homozygosity: %s (LOH-like but < %.0f%% of the median het-SNP density, not called)" % (
            ", ".join("%s:%d-%d" % (d["chrom"], d["start"], d["end"]) for d in s["germline_roh"]) or "-",
            100 * args.het_density_min),
        "subclonal_events: %s" % (("%d (>= %.0f Mb)" % (len(s["sub"]), args.subclonal_min_len / 1e6))
                                  if s["have_raw"] else "not assessed (no segments_raw)")
        + "".join("; %s:%d-%d %s ~%.0f%% of %s" % (d["chrom"], d["start"], d["end"], d["kind"],
                                                  100 * share(min(d["cn"], 1.0)), cells) for d in s["sub"]),
        "verdict: %s" % s["verdict"],
    ]
    lines += expected_lines(s, args)
    return lines


def expected_lines(s, args):
    """With --expect_gains / --expect_loh: gains found, missed and unexpected; LOH fraction in each region."""
    gains = [g.strip() for g in args.expect_gains.split(",") if g.strip()]
    if not (gains or args.expect_loh):
        return []
    lines = ["expected_gains: found %s | missed %s | unexpected %s" % (
        ",".join(g for g in gains if g in s["gained"]) or "-",
        ",".join(g for g in gains if g not in s["gained"]) or "-",
        ",".join(c for c in s["gained"] if c not in gains) or "-")]
    from .genome import parse_regions
    for chrom, start, end in parse_regions(args.expect_loh):
        tot = lohl = 0
        for d in s["segs"]:
            ov = min(d["end"], end) - max(d["start"], start) + 1
            if d["chrom"] == chrom and ov > 0:
                tot += ov
                lohl += ov if d["c_minor"] == 0 else 0
        lines.append("expected_loh %s:%d-%d: LOH fraction %.2f" % (chrom, start, end,
                                                                   lohl / tot if tot else float("nan")))
    return lines


def write_tsv(path, header, rows):
    with open(path, "w") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join(r) + "\n")


def write_call_tables(s, prefix):
    """<prefix>.chrom.tsv, .loh.tsv, .subclonal.tsv, .gains.tsv."""
    pc, pt = s["per_chrom"], s["per_chrom_tool"]
    write_tsv(prefix + ".chrom.tsv", ["chrom", "mean_cn", "loh_fraction", "modal_state", "modal_state_fraction",
                                      "ascat_modal_state"],
              [[c, "%.2f" % pc[c]["cn"], "%.3f" % pc[c]["loh"], pc[c]["state"], "%.3f" % pc[c]["state_frac"],
                pt[c]["state"]] for c in CHROMS if c in pc])
    write_tsv(prefix + ".loh.tsv", ["chrom", "start", "end", "length_mb", "mean_cn", "type", "tumour_cell_fraction"],
              [[d["chrom"], "%d" % d["start"], "%d" % d["end"], "%.2f" % (length(d) / 1e6), "%.2f" % d["cn"],
                d["type"], "%.2f" % d["frac"]] for d in s["loh"]])
    write_tsv(prefix + ".subclonal.tsv", ["chrom", "start", "end", "length_mb", "event", "background_state", "major",
                                          "minor", "major_raw", "minor_raw", "tumour_cell_fraction"],
              [[d["chrom"], "%d" % d["start"], "%d" % d["end"], "%.2f" % (length(d) / 1e6), d["kind"], d["state"]]
               + ["%.2f" % x for x in (d["a"], d["b"], d["a_raw"], d["b_raw"], min(d["cn"], 1.0))] for d in s["sub"]])
    write_tsv(prefix + ".gains.tsv", ["chrom", "start", "end", "length_mb", "state", "major", "minor",
                                      "tumour_cell_fraction"],
              [[d["chrom"], "%d" % d["start"], "%d" % d["end"], "%.2f" % (length(d) / 1e6), d["state"]]
               + ["%.2f" % x for x in (d["a"], d["b"], d["frac"])] for d in s["seg_gains"]])


def standard_output(s, args):
    """The cn_summary.json document (schema SCHEMA_VERSION)."""
    segs, fit, gate, ys = s["segs"], s["fit"], s["gate"], s["y"]
    arm_rows = arm_table(segs, args.sex, s["arms"]) if segs else []
    events = build_events(arm_rows, s["loh"], s["sub"], args.clonal_min, s["cell_share"], s["tcf"],
                          args.subclonal_max, s["seg_gains"]) if segs else []
    heh = args.heh_verdict == "yes" and segs
    n = s["count"]
    sample = dict(sample=os.path.basename(args.prefix), tool=args.label, tool_version=args.tool_version or None,
                  sex=args.sex, purity=_num4(s["purity_reported"]), purity_status=s["purity_status"],
                  purity_lower_bound=s["purity_min"], purity_fit=_num4(s["purity"]),
                  ploidy=_num4(fit.get("ploidy")), goodness_of_fit=_num4(fit.get("goodness_of_fit")),
                  aberrant_het_fraction=_num4(fit.get("aberrant_het_fraction")),
                  purity_floor=_num4(fit.get("purity_floor")),
                  chromosome_count=n if segs else None, count_class=count_class(n) if segs else None,
                  chromosome_count_range=list(s["range"]) if s["range"] else None,
                  borderline_chromosomes=sorted(s["borderline"], key=lambda c: ORDER[c]),
                  chromosome_count_uncorrected=s["count_tool"],
                  heh_key_trisomies=s["key_tri"] if heh else None, triple_trisomy=s["triple"] if heh else None,
                  sex_data=s["sex_data"] or None, sex_check=s["sex_status"] or "not_assessed",
                  x_het_ratio=None if s["xr"] is None else round(s["xr"], 3),
                  tumour_in_normal=s["tin"], tumour_in_normal_se=_num4(fit.get("tin_se")),
                  snv_purity=_num4(gate.get("snv_purity")) if gate else None,
                  snv_purity_interval=[_num4(gate.get("snv_purity_lo")), _num4(gate.get("snv_purity_hi"))]
                  if gate else None,
                  refit_gate=gate.get("decision") if gate else None,
                  y_copies=ys["y_copies"] if ys else None, y_lost_fraction=ys["lost"] if ys else None,
                  gained_chromosomes=s["gained"], lost_chromosomes=s["lost"],
                  verdict_code=s["verdict_code"] if segs else "not_assessable",
                  verdict=s["verdict"] if segs else "no segments; not assessed",
                  fractions_are="tumour cells" if s["purity_reported"] is not None else "all cells",
                  flags=qc.flags(s, fit, gate, args))
    artefacts = [dict(type="possible_germline_homozygosity", chrom=d["chrom"], start=d["start"], end=d["end"])
                 for d in s["germline_roh"]]
    artefacts += [dict(type="masked_region", chrom=c, start=st, end=e) for c, st, e in s["mask"]]
    artefacts += [dict(type="balanced_correction", chrom=d["chrom"], start=d["start"], end=d["end"],
                       state_from="%d+%d" % (d["major"], d["minor"]), state_to="%d+%d" % (d["c_major"], d["c_minor"]))
                  for d in s["corrections"]]
    arms = []
    for r in arm_rows:
        r = dict(r, **{"class": r["cls"]})
        del r["cls"]
        arms.append(r)
    return dict(schema_version=SCHEMA_VERSION, sample=sample, params=settings.params(args),
                conventions=dict(reference="GRCh38", coordinates="1-based inclusive, no chr prefix",
                                 copies="per tumour cell (unrounded *_raw) and integer state major+minor",
                                 cell_share="fraction of all cells = tumour_cell_fraction x purity",
                                 expected_copies="2; 1 for X in XY",
                                 chromosome_count=COUNT_BASIS,
                                 detection_limit="subclonal events below ~%.0f%% (loss/LOH) / %.0f%% (gain) of cells "
                                                 "are not called" % (100 * args.subclonal_min_dev_loh,
                                                                     100 * args.subclonal_min_dev)),
                arms=arms, events=events, artefacts=artefacts)


ARM_COLS = ["chrom", "arm", "start", "end", "modal_start", "modal_end", "expected", "state", "total", "arm_fraction",
            "class", "loh_fraction", "major_raw", "minor_raw", "n_het"]
EVENT_COLS = ["scope", "type", "chrom", "start", "end", "length_mb", "arms", "state", "copies", "major_raw",
              "minor_raw", "clonality", "tumour_cell_fraction", "cell_share"]


def report_text(lines):
    """The report as printed and written to heh.txt: each part of a long value on its own indented line."""
    return "\n".join(lines).replace("; ", "\n    ") + "\n"


def write_all(s, args):
    """Every output file; returns the report text."""
    prefix = "%s.%s" % (args.prefix, args.label)
    write_call_tables(s, prefix)
    text = report_text(text_lines(s, args))
    with open(prefix + ".heh.txt", "w") as f:
        f.write(text)
    doc = standard_output(s, args)
    with open(prefix + ".cn_summary.json", "w") as f:
        json.dump(doc, f, indent=1)

    def cell(v):
        return "NA" if v is None else str(v)
    write_tsv(prefix + ".arms.tsv", ARM_COLS, [[cell(r[k]) for k in ARM_COLS] for r in doc["arms"]])
    write_tsv(prefix + ".events.tsv", EVENT_COLS, [[cell(e[k]) for k in EVENT_COLS] for e in doc["events"]])
    return text
