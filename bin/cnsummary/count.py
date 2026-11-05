"""The chromosome count, karyotype style, its range, and the HeH verdict.

A chromosome counts as gained or lost only when both arms agree (acrocentric: the q arm); arm-level changes are
listed, not counted. The count is copy-number based, not a karyotype (COUNT_BASIS).
"""
from . import settings
from .genome import ACROCENTRIC, AUTOSOMES, expected_cn, length

COUNT_BASIS = ("copy-number based: the modal copies of each chromosome summed; not a karyotype. Structural, "
               "marker and ring chromosomes, Robertsonian translocations and acrocentric p arms are not seen, "
               "clones are averaged, and Y is assumed present for XY")


def chromosome_copies(arms, sex, whole_frac):
    """Copies of each chromosome (autosomes and X): the modal total when both arms share it over >= half of each
    arm (acrocentric: the q arm over >= whole_frac), else the expected copies."""
    out = {}
    for c in AUTOSOMES + ["X"]:
        p, q = arms.get((c, "p")), arms.get((c, "q"))
        exp = expected_cn(c, sex)
        if c in ACROCENTRIC or p is None:
            out[c] = q[0] if q and q[1] >= whole_frac else exp
        elif q is None:
            out[c] = exp
        else:
            out[c] = p[0] if p[0] == q[0] and p[1] >= settings.ARM_MIN and q[1] >= settings.ARM_MIN else exp
    return out


def borderline_chromosomes(arms, copies, sex, whole_frac, margin=0.05):
    """Chromosomes whose count decision is within +-margin of its coverage threshold (whole_frac for an acrocentric
    q arm, half of each arm otherwise): {chrom: (counted, alternative)}."""
    out = {}
    for c in AUTOSOMES + ["X"]:
        p, q = arms.get((c, "p")), arms.get((c, "q"))
        exp = expected_cn(c, sex)
        if q is None or q[0] == exp:
            continue
        if c in ACROCENTRIC or p is None:
            near = abs(q[1] - whole_frac) < margin
        else:
            near = p[0] == q[0] and abs(min(p[1], q[1]) - settings.ARM_MIN) < margin
        if near:
            out[c] = (copies[c], exp if copies[c] != exp else q[0])
    return out


def arm_changes(arms, copies):
    """Arms whose modal total (over >= half the arm) differs from their chromosome's counted copies, e.g.
    "7q loss (1+0)"."""
    out = []
    for c in AUTOSOMES + ["X"]:
        for arm in ("p", "q"):
            v = arms.get((c, arm))
            if v is None or (arm == "p" and c in ACROCENTRIC) or v[1] < settings.ARM_MIN:
                continue
            if v[0] != copies[c]:
                out.append("%s%s %s (%s)" % (c, arm, "gain" if v[0] > copies[c] else "loss", v[2]))
    return out


def chromosome_count(copies, sex, y_present=True):
    """Copies of the autosomes and X, plus Y for a male unless Y is lost in most tumour cells."""
    return sum(copies.values()) + (1 if sex == "XY" and y_present else 0)


def purity_dependent_gains(copies, sub, segs, sex, whole_frac):
    """When purity is only bounded or unknown, a gain that is clonal at a low purity is subclonal at purity 1 and
    the other way round, so whether it counts cannot be decided. -> [(chrom, counted)]: the chromosomes counted
    as gained, and those whose subclonal gain covers >= whole_frac of their length."""
    out = []
    for c in AUTOSOMES + ["X"]:
        exp = expected_cn(c, sex)
        if copies.get(c, exp) > exp:
            out.append((c, True))
            continue
        total = sum(length(d) for d in segs if d["chrom"] == c)
        gain = sum(length(d) for d in sub if d["chrom"] == c and d["kind"] == "subclonal gain")
        if total and gain >= whole_frac * total:
            out.append((c, False))
    return out


def count_range(n, n_tool, borderline, purity_alt, copies, sex, y_present, y_border):
    """(low, high) around the count n: borderline chromosomes either way, the tool's own count before the
    balance correction, gains whose counting depends on purity, and Y when lost in about half the tumour cells."""
    lo = n + sum(min(alt - cnt, 0) for cnt, alt in borderline.values())
    hi = n + sum(max(alt - cnt, 0) for cnt, alt in borderline.values())
    lo, hi = min(lo, n_tool), max(hi, n_tool)
    lo = min(lo, n - sum(copies[c] - expected_cn(c, sex) for c, counted in purity_alt if counted))
    hi = max(hi, n + sum(1 for c, counted in purity_alt if not counted))
    if y_border:
        lo, hi = lo - (1 if y_present else 0), hi + (0 if y_present else 1)
    return lo, hi


def heh_trisomies(copies):
    """{chrom: gained} for the prognostic HeH chromosomes, and whether +4, +10 and +17 are all present."""
    key = {c: copies.get(c, 2) >= 3 for c in settings.HEH_KEY_TRISOMIES}
    return key, all(key[c] for c in ("4", "10", "17"))


def verdict(n, whole_loh, heh_verdict, masked_hypo_min_loh):
    """(code, text). Doubled (masked) hypodiploidy can reach an HeH count; it shows as many whole autosomes in
    LOH (2+0), so it is flagged before the count is read."""
    if heh_verdict == "no":
        return "not_applicable", "not applicable (HeH verdict is for ALL only)"
    if len(whole_loh) >= masked_hypo_min_loh:
        return ("possible_doubled_hypodiploidy",
                "POSSIBLE DOUBLED (MASKED) HYPODIPLOIDY — %d whole-autosome LOH >= %d; review before reporting HeH"
                % (len(whole_loh), masked_hypo_min_loh))
    if settings.HEH_RANGE[0] <= n <= settings.HEH_RANGE[1]:
        return "HeH", "consistent with HeH"
    return "not_HeH", "not HeH by chromosome count"


def count_class(n):
    """Modal-number classes as used for ALL (near-haploid and low hypodiploid are the high-risk groups)."""
    for hi, name in ((30, "near-haploid (24-30)"), (39, "low hypodiploid (31-39)"), (43, "high hypodiploid (40-43)"),
                     (45, "near-diploid hypodiploid (44-45)"), (46, "diploid (46)"), (50, "low hyperdiploid (47-50)"),
                     (67, "high hyperdiploid (51-67)"), (80, "near-triploid (68-80)")):
        if n <= hi:
            return name
    return "near-tetraploid (>80)"
