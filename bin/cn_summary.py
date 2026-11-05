#!/usr/bin/env python3
"""Copy-number summary of one tumour/normal fit: LOH, gains, subclonal events, chromosome count, HeH verdict
(stdlib only). The same summary for every tool: ASCAT's tables directly, FACETS' after to_segments.py.

Inputs  <fit>.purity_ploidy.tsv   purity, ploidy, purity_status, QC columns (sex, tumour in normal, purity floor)
        <fit>.segments.tsv        chr startpos endpos nMajor nMinor                  (integer states)
        <fit>.segments_raw.tsv    ... nAraw nBraw n_het                              (unrounded copies; optional)
Outputs <prefix>.<label>.heh.txt          the report, one "key: value" line per item
        <prefix>.<label>.cn_summary.json  the standardised summary (sample, arms, events, artefacts, params)
        <prefix>.<label>.{arms,events}.tsv          the arm and event layers of the JSON
        <prefix>.<label>.{chrom,loh,gains,subclonal}.tsv   per chromosome, LOH, segmental gains, subclonal events

How the calls are made (details in the docstrings of cnsummary/):
  - Unrounded copies are used, as integer states hide subclonal events and can give a balanced segment an odd
    state. A purity under-fit inflates every allele's distance from one copy by the same factor; it is estimated
    from the 2+1 segments and divided out. Balanced segments get an even state.
  - The chromosome count is karyotype style: a chromosome is gained or lost only when both arms agree
    (acrocentric: the q arm); arm-level changes are listed, not counted. Arms are split at the centromere region.
    A borderline chromosome, the tool's own count, a purity-dependent gain or a half-lost Y make it a range.
  - Artefact filters: a segment with far fewer het SNPs than the median is germline homozygosity, not LOH;
    --mask regions (default the X-transposed region) are left out of LOH and subclonal calls.
  - When purity is not identifiable (a flat genome), fractions are reported as shares of all cells.
  - HeH: 51-67 chromosomes; +4, +10, +17, +18 are reported. Doubled hypodiploidy can reach the same count and
    shows as many whole autosomes in LOH (2+0); it is flagged first.
All thresholds (cnsummary/settings.py) are heuristics, not validated.
"""
import argparse

from cnsummary import settings
from cnsummary.report import write_all
from cnsummary.summary import summarise


def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--purity_ploidy", required=True)
    ap.add_argument("--segments", required=True)
    ap.add_argument("--segments_raw", default="", help="unrounded copies (enables the correction and subclonal calls)")
    ap.add_argument("--prefix", required=True)
    ap.add_argument("--label", default="ascat", help="the tool; outputs are <prefix>.<label>.*")
    ap.add_argument("--sex", choices=["XX", "XY"], required=True)
    ap.add_argument("--ydepth", default="", help="Y and autosome read counts, tumour and normal (y_depth.py); XY only")
    ap.add_argument("--gate", default="", help="refit decision from snv_purity.py --gate_out (reported and flagged)")
    ap.add_argument("--tool_version", default="", help="recorded in the JSON")
    ap.add_argument("--expect_gains", default="", help="known gains, comma-separated chromosomes (checked)")
    ap.add_argument("--expect_loh", default="", help="known LOH, chrom:start-end[,...] (checked)")
    settings.add_arguments(ap)
    return ap.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    print(write_all(summarise(args), args), end="")


if __name__ == "__main__":
    main()
