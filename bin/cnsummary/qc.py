"""Quality checks of the pair and the fit: purity status, sex, tumour in the normal, loss of Y, and the flags."""
import os

from . import settings
from .genome import AUTOSOMES, CENTROMERE, SEX_XY_MAX, expected_cn, length, y_autosome_windows
from .tables import num, read_tsv


def purity_state(status, pur, segs, loh, sub, sex, min_len):
    """-> (purity_status, lower_bound). status comes from the fit (ok / not identifiable / manual).
    lower_bound: there are gains but no loss or LOH anywhere, so a clonal gain at purity p fits as well as a
    subclonal gain in a share p of the cells at purity 1; purity is then only >= the largest gained share."""
    status = {"not identifiable": "not_identifiable"}.get(status, status)
    if not segs:
        return "no_fit", None
    if pur is None:
        return ("not_identifiable" if status == "ok" else status), None
    if status != "ok":
        return status, None
    big = [d for d in segs if length(d) >= min_len and not d.get("masked") and not d.get("lowhet")
           and d["chrom"] in CENTROMERE]
    lossy = loh or any(s["kind"] != "subclonal gain" for s in sub) or \
        any(d["c_major"] + d["c_minor"] < expected_cn(d["chrom"], sex) for d in big)
    shares = [min(s["cn"], 1.0) * pur for s in sub if s["kind"] == "subclonal gain"] + \
             [pur for d in big if d["c_major"] + d["c_minor"] > expected_cn(d["chrom"], sex)]
    if lossy or not shares:
        return status, None
    return "lower_bound", round(max(shares), 3)


def x_het_ratio(segs):
    """Het-SNP density on X relative to the autosomes (het SNPs per Mb). The het SNPs are germline, so tumour
    changes do not move it: ~0.6-1 with two X copies, ~0 with one. None without n_het or with < 50 Mb."""
    def dens(keep):
        ds = [d for d in segs if keep(d) and d["n_het"] is not None and not d.get("masked")]
        n = sum(length(d) for d in ds)
        return (sum(d["n_het"] for d in ds) / (n / 1e6)) if n >= 50e6 else None
    x, auto = dens(lambda d: d["chrom"] == "X"), dens(lambda d: d["chrom"] in AUTOSOMES)
    return None if x is None or not auto else x / auto


def sex_check(fit, sex, xr):
    """(sex from the data, status ok / mismatch / unclear / ""). The fit's check (the normal's X heterozygosity)
    is used when it was made; the X het-SNP density of the segments (xr) adds a check that works for any tool and
    catches a male given as XX."""
    def blank(v):
        return "" if v in (None, "", "NA", "not_assessed", "unknown") else v
    data, status = blank(fit.get("sex_data")), blank(fit.get("sex_check"))
    if not status and xr is not None and sex == "XX":
        data, status = ("XY", "mismatch") if xr < SEX_XY_MAX else ("XX", "ok")
    elif status == "ok" and sex == "XX" and xr is not None and xr < SEX_XY_MAX:
        status = "mismatch"
    return data, status


def y_status(path, segs, pur, min_loss=settings.Y_LOSS_MIN):
    """Loss of Y from read counts (y_depth.py): r = (Y / autosome windows) in the tumour over the normal.
    Normal cells carry one Y and two autosome copies; tumour cells y Y copies and the windows' copies c (from the
    unrounded segments). Tumour sample: Y depth ~ p*y + (1 - p), window depth ~ p*c + 2(1 - p); normal: 1 and 2, so
    y = (r * mean(p*c + 2(1 - p)) / 2 - (1 - p)) / p. Without a purity: share of all cells without Y ~ 1 - r.
    -> dict(ratio, y_copies, lost, cells, called) or None without counts."""
    if not path or not os.path.exists(path):
        return None
    rows = {r["role"]: r for r in read_tsv(path)}
    try:
        (yt, at), (yn, an) = [(float(rows[k]["y_reads"]), float(rows[k]["auto_reads"])) for k in ("tumor", "normal")]
        r = (yt / at) / (yn / an)
    except (KeyError, ValueError, ZeroDivisionError):
        return None

    def copies(c, pos):
        return next((d["a"] + d["b"] for d in segs if d["chrom"] == c and d["start"] <= pos <= d["end"]), 2.0)
    wins = y_autosome_windows()
    if pur:
        m = sum(pur * copies(c, (s + e) // 2) + 2 * (1 - pur) for c, s, e in wins) / len(wins)
        y = (r * m / 2 - (1 - pur)) / pur
        lost, cells = 1 - min(max(y, 0.0), 1.0), "tumour cells"
    else:
        y, lost, cells = None, 1 - min(max(r, 0.0), 1.0), "all cells"
    return dict(ratio=round(r, 3), y_copies=None if y is None else round(y, 2), lost=round(lost, 3),
                cells=cells, called=lost >= min_loss)


def flags(s, fit, gate, args):
    """The flags of cn_summary.json, from the computed summary s."""
    out = []
    pstatus, pur, pur_rep = s["purity_status"], s["purity"], s["purity_reported"]
    if pstatus == "not_identifiable":
        out.append("purity_not_identifiable")
    if pstatus == "lower_bound":
        out.append("purity_lower_bound_only")
    if pur_rep is not None and pur_rep < args.low_purity:
        out.append("low_tumour_content")
    if fit.get("fit_used") == "allele":
        out.append("purity_from_both_allele_refit")
    if s["scale"] > settings.SCALE_FLAG:
        out.append("gained_alleles_rescaled")
    if not s["have_nhet"]:
        out.append("germline_homozygosity_filter_off")
    if args.sex == "XY" and s["y"] is None:
        out.append("y_not_assessed")
    if s["y"] and s["y"]["called"]:
        out.append("loss_of_y")
    if s["borderline"]:
        out.append("count_borderline")
    if s["count"] is not None and s["count"] != s["count_tool"]:
        out.append("count_changed_by_balance_correction")
    if s["sex_status"] == "mismatch":
        out.append("sex_mismatch")
    if s["purity_alt"]:
        out.append("count_depends_on_purity")
    if s["floor"] is not None and pur is not None and pur < s["floor"] - settings.FLOOR_MARGIN:
        out.append("purity_below_floor")
    if gate.get("decision") == "reject_refit":
        out.append("refit_rejected_by_snv")
    elif gate and gate.get("decision") not in ("accept_refit", "no_refit"):
        out.append("refit_unconfirmed")
    if s["tin"] is not None and s["tin"] >= args.tin_max:
        out.append("tumour_in_normal")
    ploidy = num(fit.get("ploidy"))
    if ploidy is not None and ploidy > settings.HIGH_PLOIDY:
        out.append("high_ploidy_check_doubling")
    return out
