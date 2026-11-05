"""GRCh38 geometry shared by every script: chromosome names and order, centromeres, arms, PARs, lengths.

Coordinates are 1-based and inclusive, contig names without "chr", as in our reference.
"""

AUTOSOMES = [str(i) for i in range(1, 23)]
CHROMS = AUTOSOMES + ["X", "Y"]
FIT_CHROMS = AUTOSOMES + ["X"]          # what ASCAT and FACETS are run on
ORDER = {c: i for i, c in enumerate(CHROMS)}

# Centromere regions in Mb: UCSC cytoBand "acen" bands, p11.1 start to q11 end. The p arm ends where the region
# starts and the q arm begins where it ends; bases inside belong to neither arm (few SNPs, unreliable calls).
CENTROMERE = {"1": (121.7, 125.1), "2": (91.8, 96.0), "3": (87.8, 94.0), "4": (48.2, 51.8), "5": (46.1, 51.4),
              "6": (58.5, 62.6), "7": (58.1, 62.1), "8": (43.2, 47.2), "9": (42.2, 45.5), "10": (38.0, 41.6),
              "11": (51.0, 55.8), "12": (33.2, 37.8), "13": (16.5, 18.9), "14": (16.1, 18.2), "15": (17.5, 20.5),
              "16": (35.3, 38.4), "17": (22.7, 27.4), "18": (15.4, 21.5), "19": (24.2, 28.1), "20": (25.7, 30.4),
              "21": (10.9, 13.0), "22": (13.7, 17.4), "X": (58.1, 63.8)}
ACROCENTRIC = {"13", "14", "15", "21", "22"}   # the p arm has no SNPs: the chromosome is its q arm

# Autosome lengths in Mb
CHROM_LEN = {"1": 248.96, "2": 242.19, "3": 198.30, "4": 190.21, "5": 181.54, "6": 170.81, "7": 159.35, "8": 145.14,
             "9": 138.39, "10": 133.80, "11": 135.09, "12": 133.28, "13": 114.36, "14": 107.04, "15": 101.99,
             "16": 90.34, "17": 83.26, "18": 80.37, "19": 58.62, "20": 64.44, "21": 46.71, "22": 50.82}

# X pseudoautosomal regions (PAR1, PAR2); X outside them is single-copy in a male
PAR_X = ((0, 2_781_479), (155_701_383, 156_030_895))

# Sex from germline heterozygosity on X (non-PAR) relative to the autosomes: two X copies give roughly the
# autosomal rate, one X gives ~0. In between: unclear.
SEX_XX_MIN, SEX_XY_MAX = 0.3, 0.1

# Male-specific Y, X-degenerate part: after the X-transposed region (~2.9-6.1 Mb), before the heterochromatin.
# Reads with MAPQ >= 30 there are Y-specific; a female sample has about none.
Y_REGION = ("Y", 6_500_000, 26_000_000)


def norm_chrom(c):
    """Contig name without "chr"; FACETS' 23 is X."""
    c = str(c)
    c = c[3:] if c.startswith("chr") else c
    return "X" if c == "23" else c


def arm_part(c, arm, start, end):
    """(start, end) of [start, end] inside the given arm of chromosome c (end <= start: none)."""
    cs, ce = (int(x * 1e6) for x in CENTROMERE[c])
    return (start, min(end, cs)) if arm == "p" else (max(start, ce), end)


def y_autosome_windows(half=500_000):
    """1 Mb windows at the middle of every autosome arm (acrocentrics: q only): the depth reference for Y."""
    out = []
    for c in AUTOSOMES:
        cs, ce = (x * 1e6 for x in CENTROMERE[c])
        mids = ([] if c in ACROCENTRIC else [cs / 2]) + [(ce + CHROM_LEN[c] * 1e6) / 2]
        out += [(c, int(m - half), int(m + half)) for m in mids]
    return out


def expected_cn(c, sex):
    """Copies of a chromosome in a normal cell: 2, or 1 for X in a male."""
    return 1 if (c == "X" and sex == "XY") else 2


def length(d):
    """Length in bp of anything with "start" and "end" (inclusive)."""
    return d["end"] - d["start"] + 1


def overlap(a, b):
    """Bases shared by two intervals on the same chromosome (0 when they do not overlap or the chromosomes
    differ). Both are dicts with "chrom", "start", "end"."""
    if a["chrom"] != b["chrom"]:
        return 0
    return max(0, min(a["end"], b["end"]) - max(a["start"], b["start"]) + 1)


def parse_regions(spec):
    """"chrom:start-end[,chrom:start-end]" -> [(chrom, start, end)]."""
    out = []
    for reg in filter(None, spec.split(",")):
        c, rng = reg.split(":")
        s, e = (int(float(x)) for x in rng.split("-"))
        out.append((norm_chrom(c), s, e))
    return out
