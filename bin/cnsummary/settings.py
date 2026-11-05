"""Every threshold of the summary, defined once: the command line and the JSON "params" are built from here.

All values are heuristics tuned on the development pairs; none is validated. config.env sets the ones the summary job passes.
"""

# (name, type, default, meaning). In this order in the JSON "params".
THRESHOLDS = [
    ("whole_frac", float, 0.9,
     "a state must cover this fraction of a chromosome (acrocentric: of its q arm) to count as whole-chromosome"),
    ("masked_hypo_min_loh", int, 5,
     "flag possible doubled hypodiploidy at >= this many whole autosomes in LOH"),
    ("min_loh_len", int, 1000000, "shortest LOH segment reported (bp)"),
    ("max_gap", int, 1000000,
     "neighbouring pieces of one event closer than this merge; LOH pieces shorter than this are dropped first (bp)"),
    ("balanced_tol", float, 0.15, "|major - minor| below this: allelically balanced (unrounded copies)"),
    ("subclonal_min_dev", float, 0.15, "a gained allele at least this far from its integer state is subclonal"),
    ("subclonal_min_dev_loh", float, 0.12,
     "a minor allele at least this far below 1 is a subclonal LOH/cnLOH (allele-ratio based, more precise)"),
    ("subclonal_min_len", int, 5000000, "shortest subclonal event or segmental gain reported (bp)"),
    ("balanced_gain_min_len", int, 20000000,
     "a subclonal gain without allelic shift (k+k: read depth only) needs this much such gain on its chromosome (bp)"),
    ("het_density_min", float, 0.25,
     "LOH needs >= this x the median het-SNP density; below it the normal is homozygous too (germline)"),
    ("mask", str, "X:88000000-93500000",
     "chrom:start-end[,...] left out of LOH and subclonal calls (default: the X-transposed region)"),
    ("low_purity", float, 0.2, "flag low_tumour_content below this purity"),
    ("clonal_min", float, 0.85, "events in >= this fraction of tumour cells are clonal"),
    ("subclonal_max", float, 0.7,
     "events in < this fraction of tumour cells are subclonal; in between: possibly_subclonal"),
    ("borderline_margin", float, 0.05,
     "a chromosome whose arm coverage is this close to its count threshold makes the count a range"),
    ("heh_verdict", str, "yes", "no: not an ALL sample; the HeH verdict is not_applicable"),
    ("tin_max", float, 0.05, "flag tumour_in_normal at >= this fraction of tumour cells in the normal"),
]

# Fixed rules that are not options
SCALE_MIN_SEGMENTS = 3          # 2+1 segments needed to estimate the allele scale
SCALE_MAX = 1.6                 # the scale is capped here
SCALE_FLAG = 1.1                # report and flag a scale above this
HIGH_PLOIDY = 3.2               # ploidy above this: check for a doubled solution
FLOOR_MARGIN = 0.05             # flag a fit purity this far below the purity floor
Y_LOSS_MIN = 0.15               # Y lost in >= this share of cells is reported
COVERED = 0.9                   # a segment event >= this much inside an arm/chromosome event of the same type is not repeated
ARM_MIN = 0.5                   # an arm's modal state must cover this much of the arm
HEH_RANGE = (51, 67)            # high hyperdiploid chromosome count
HEH_KEY_TRISOMIES = ("4", "10", "17", "18")   # prognostic in HeH ALL; +4 +10 +17 together: triple trisomy


def add_arguments(ap):
    """Adds every threshold to an argparse parser as --<name>."""
    for name, typ, default, meaning in THRESHOLDS:
        kw = dict(type=typ, default=default, help=meaning.replace("%", "%%"))
        if name == "heh_verdict":
            kw["choices"] = ["yes", "no"]
        ap.add_argument("--" + name, **kw)


def params(args):
    """The thresholds used, for the JSON."""
    return {name: getattr(args, name) for name, _, _, _ in THRESHOLDS}
