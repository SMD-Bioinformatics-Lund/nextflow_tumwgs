#!/usr/bin/env python
"""Convert Pelops JSON output to a Manta-compatible BND VCF.

Usage:
    pelops_to_vcf.py \\
        --tumor-json  tumor_pelops.json \\
        --normal-json normal_pelops.json \\
        --tumor-name  TUMOR_SAMPLE \\
        --normal-name NORMAL_SAMPLE \\
        --output      pelops.vcf

Filter logic:
  - CoreDUX4-IGH (id=01) and ExtendedDUX4-IGH (id=02):
        PASS if T_SRPB >= srpb_threshold)
        LowSRPB otherwise
    Fold-change is NOT applied to IGH-DUX4 pairs. Normal blood contains
    circulating B lymphocytes undergoing V(D)J recombination at the IGH
    locus, which generates inherent DUX4-IGH discordant read background.
    Absolute SRPB is the validated threshold for these calls (Pelops paper).

  - UNNAMED rearrangements (id >= 03):
        PASS if T_SRPB >= srpb_threshold AND ln(T_SRPB / N_SRPB) > log_fold_threshold
        NormalSignal if T_SRPB >= srpb_threshold but fold-change insufficient
        LowSRPB if T_SRPB < srpb_threshold

  Use --strict-fold-change to apply fold-change to all rearrangement types.

Coordinate strategy (Option A — simple midpoints):
  POS is the midpoint of the representative sub-region. For DUX4 regions the
  most-telomeric chr4 sub-region is used; for IGH the single region midpoint.
  Positional uncertainty is encoded in CIPOS/CIEND as ±half-width of the
  representative region. POS values are NOT assembly-resolved breakpoints.

FORMAT fields (Manta convention — no GT field):
  PR : ref_count,alt_count  — ref is set to 0 (not computed by Pelops)
  SR : ref_count,alt_count  — ref is set to 0 (not computed by Pelops)
  FORMAT is PR:SR when split reads > 0, else PR only (matches Manta behaviour)

Chromosome naming:
  Output uses no-chr contig names (1, 2, ..., X, Y, M) to match Manta's
  GCA_000001405.15_GRCh38_no_alt_analysis_set_nochr reference.
"""

import argparse
import json
import math
import sys
from typing import Dict, List, Optional, Tuple

# GRCh38 primary assembly contig lengths — no-chr names (matches Manta nochr reference)
GRCH38_CONTIGS: Dict[str, int] = {
    "1": 248956422,
    "2": 242193529,
    "3": 198295559,
    "4": 190214555,
    "5": 181538259,
    "6": 170805979,
    "7": 159345973,
    "8": 145138636,
    "9": 138394717,
    "10": 133797422,
    "11": 135086622,
    "12": 133275309,
    "13": 114364328,
    "14": 107043718,
    "15": 101991189,
    "16": 90338345,
    "17": 83257441,
    "18": 80373285,
    "19": 58617616,
    "20": 64444167,
    "21": 46709983,
    "22": 50818468,
    "X": 156040895,
    "Y": 57227415,
    "M": 16569,
}

# Chromosome sort order for deterministic output
_CHROM_ORDER = {c: i for i, c in enumerate(GRCH38_CONTIGS)}


def _strip_chr(chrom: str) -> str:
    """Strip 'chr' prefix to match no-chr Manta reference (chrM → M, chrX → X, etc.)."""
    return chrom[3:] if chrom.startswith("chr") else chrom


def _chrom_key(chrom: str) -> int:
    return _CHROM_ORDER.get(_strip_chr(chrom), 999)


def pick_representative_region(regions: List[dict]) -> dict:
    """Choose a single sub-region to represent a CompoundRegion.

    Priority:
      1. chr4 sub-regions, sorted by start descending (most-telomeric DUX4)
      2. chr10 sub-regions, sorted by start descending
      3. Any chromosome: largest region by span
    """
    chr4 = sorted(
        [r for r in regions if _strip_chr(r["chrom"]) == "4"],
        key=lambda r: r["start"],
        reverse=True,
    )
    if chr4:
        return chr4[0]

    chr10 = sorted(
        [r for r in regions if _strip_chr(r["chrom"]) == "10"],
        key=lambda r: r["start"],
        reverse=True,
    )
    if chr10:
        return chr10[0]

    return max(regions, key=lambda r: r["end"] - r["start"])


def midpoint(region: dict) -> int:
    return (region["start"] + region["end"]) // 2


def halfwidth(region: dict) -> int:
    return (region["end"] - region["start"]) // 2


def regions_to_str(regions: List[dict]) -> str:
    """Encode sub-regions as pipe-delimited string for INFO field (no-chr names)."""
    sorted_regions = sorted(regions, key=lambda r: (_chrom_key(r["chrom"]), r["start"]))
    return "|".join(f"{_strip_chr(r['chrom'])}:{r['start']}-{r['end']}" for r in sorted_regions)


def is_igh_dux4(a_name: str, b_name: str) -> bool:
    """True for the two named IGH-DUX4 rearrangement types (id 01 and 02)."""
    return b_name == "IGH" and a_name in ("CoreDUX4", "ExtendedDUX4")


def compute_filter(
    a_name: str,
    b_name: str,
    t_srpb: float,
    n_srpb: float,
    srpb_threshold: float,
    log_fold_threshold: float,
    strict_fold_change: bool,
) -> str:
    if t_srpb < srpb_threshold:
        return "LowSRPB"

    # For IGH-DUX4 pairs: fold-change not applied unless --strict-fold-change
    if is_igh_dux4(a_name, b_name) and not strict_fold_change:
        return "PASS"

    # UNNAMED or strict mode: apply fold-change filter
    if n_srpb == 0.0:
        return "PASS"

    fold = math.log(t_srpb / n_srpb)  # natural log
    return "PASS" if fold > log_fold_threshold else "NormalSignal"


def bnd_alts(chrom_a: str, pos_a: int, chrom_b: str, pos_b: int) -> Tuple[str, str]:
    """Return (alt_for_record_A, alt_for_record_B) in BND notation.

    Orientation convention: forward-forward translocation.
      Record A ALT: N[B:posB[
      Record B ALT: ]A:posA]N
    Contig names are no-chr (matching Manta reference).
    """
    alt_a = f"N[{chrom_b}:{pos_b}["
    alt_b = f"]{chrom_a}:{pos_a}]N"
    return alt_a, alt_b


def make_bnd_pair(
    r_tumor: dict,
    r_normal: Optional[dict],
    normal_name: str,
    tumor_name: str,
    srpb_threshold: float,
    log_fold_threshold: float,
    strict_fold_change: bool,
) -> List[str]:
    """Generate a VCF BND record pair for one rearrangement."""
    ev_t = r_tumor["evidence"]
    ev_n = r_normal["evidence"] if r_normal else {"paired_reads": 0, "split_reads": 0, "SRPB": 0.0}

    t_srpb = ev_t["SRPB"]
    n_srpb = ev_n["SRPB"]
    a_name = r_tumor["A"]["name"]
    b_name = r_tumor["B"]["name"]
    rid = r_tumor["id"]

    vcf_filter = compute_filter(
        a_name, b_name, t_srpb, n_srpb, srpb_threshold, log_fold_threshold, strict_fold_change
    )

    regions_a = r_tumor["A"]["regions"]
    regions_b = r_tumor["B"]["regions"]

    rep_a = pick_representative_region(regions_a)
    rep_b = pick_representative_region(regions_b)

    chrom_a = _strip_chr(rep_a["chrom"])
    chrom_b = _strip_chr(rep_b["chrom"])
    pos_a, pos_b = midpoint(rep_a), midpoint(rep_b)
    hw_a, hw_b = halfwidth(rep_a), halfwidth(rep_b)

    id0 = f"PelopsCall:{rid}:0"
    id1 = f"PelopsCall:{rid}:1"

    alt_a, alt_b = bnd_alts(chrom_a, pos_a, chrom_b, pos_b)

    regions_a_str = regions_to_str(regions_a)
    regions_b_str = regions_to_str(regions_b)

    def make_info(mateid: str, cipos_hw: int, ciend_hw: int) -> str:
        return (
            f"SVTYPE=BND"
            f";MATEID={mateid}"
            f";CIPOS=-{cipos_hw},{cipos_hw}"
            f";CIEND=-{ciend_hw},{ciend_hw}"
            f";SOMATIC"
            f";SVMETHOD=pelops"
            f";SRPB_T={t_srpb}"
            f";SRPB_N={n_srpb}"
            f";DUX4_REGIONS_A={regions_a_str}"
            f";DUX4_REGIONS_B={regions_b_str}"
        )

    info0 = make_info(id1, hw_a, hw_b)
    info1 = make_info(id0, hw_b, hw_a)

    # Use PR:SR when tumor has split reads, PR only otherwise (matches Manta behaviour)
    has_split = ev_t["split_reads"] > 0
    if has_split:
        fmt = "PR:SR"
        normal_fmt = f"0,{ev_n['paired_reads']}:0,{ev_n['split_reads']}"
        tumor_fmt  = f"0,{ev_t['paired_reads']}:0,{ev_t['split_reads']}"
    else:
        fmt = "PR"
        normal_fmt = f"0,{ev_n['paired_reads']}"
        tumor_fmt  = f"0,{ev_t['paired_reads']}"

    def record(chrom, pos, rec_id, alt, info) -> str:
        return "\t".join([
            chrom, str(pos), rec_id, "N", alt, ".", vcf_filter,
            info, fmt, normal_fmt, tumor_fmt,
        ])

    return [
        record(chrom_a, pos_a, id0, alt_a, info0),
        record(chrom_b, pos_b, id1, alt_b, info1),
    ]


def build_header(
    normal_name: str,
    tumor_name: str,
    version: str,
    cli_command: str,
    srpb_threshold: float,
    log_fold_threshold: float,
    strict_fold_change: bool,
) -> str:
    contig_lines = "\n".join(
        f"##contig=<ID={c},length={length}>" for c, length in GRCH38_CONTIGS.items()
    )
    fold_note = (
        "fold-change applied to ALL rearrangement types (--strict-fold-change)"
        if strict_fold_change
        else "fold-change applied to UNNAMED rearrangements only; "
             "IGH-DUX4 calls use absolute SRPB threshold (V(D)J background in normal blood)"
    )
    return (
        "##fileformat=VCFv4.1\n"
        f"##source=pelops {version}\n"
        "##reference=GRCh38\n"
        "##NOTE=POS values are region midpoints not assembly-resolved breakpoints; "
        "positional uncertainty is encoded in CIPOS/CIEND\n"
        "##NOTE=PR and SR ref counts are 0 (reference-supporting reads not computed by Pelops)\n"
        f"##NOTE=Filter logic: SRPB_threshold={srpb_threshold}; "
        f"ln_fold_threshold={log_fold_threshold}; {fold_note}\n"
        f"##NOTE=pelops_cmd={cli_command}\n"
        f"{contig_lines}\n"
        '##FILTER=<ID=PASS,Description="All filters passed">\n'
        f'##FILTER=<ID=LowSRPB,Description="Tumor SRPB below threshold ({srpb_threshold})">\n'
        f'##FILTER=<ID=NormalSignal,Description="Tumor SRPB meets threshold but '
        f'ln(T_SRPB/N_SRPB) <= {log_fold_threshold} (applied to UNNAMED rearrangements)">\n'
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n'
        '##INFO=<ID=MATEID,Number=1,Type=String,Description="ID of the mate breakend record">\n'
        '##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS '
        '(half-width of representative DUX4 sub-region)">\n'
        '##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around mate POS '
        '(half-width of partner region)">\n'
        '##INFO=<ID=SVMETHOD,Number=1,Type=String,Description="Caller used to identify this variant">\n'
        '##INFO=<ID=SRPB_T,Number=1,Type=Float,Description="Tumor spanning read pairs per billion (Pelops)">\n'
        '##INFO=<ID=SRPB_N,Number=1,Type=Float,Description="Normal spanning read pairs per billion (Pelops)">\n'
        '##INFO=<ID=DUX4_REGIONS_A,Number=1,Type=String,Description="DUX4 CompoundRegion A sub-regions '
        'used by Pelops, pipe-delimited chrom:start-end">\n'
        '##INFO=<ID=DUX4_REGIONS_B,Number=1,Type=String,Description="Partner CompoundRegion B sub-regions '
        'used by Pelops, pipe-delimited chrom:start-end">\n'
        '##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="Indicates a somatic structural variant">\n'
        '##FORMAT=<ID=PR,Number=2,Type=Integer,Description="Spanning paired-end read support: '
        'ref_count,alt_count">\n'
        '##FORMAT=<ID=SR,Number=2,Type=Integer,Description="Split-read support: ref_count,alt_count">\n'
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{normal_name}\t{tumor_name}"
    )


def load_normal_index(normal_data: dict) -> Dict[Tuple[str, str], dict]:
    """Index normal rearrangements by (A_name, B_name) for O(1) lookup."""
    index = {}
    for r in normal_data.get("rearrangements", []):
        key = (r["A"]["name"], r["B"]["name"])
        index[key] = r
    return index


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Convert Pelops JSON output to Manta-compatible BND VCF",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--tumor-json",  required=True, metavar="FILE",
                        help="Pelops results JSON from the tumor sample")
    parser.add_argument("--normal-json", required=True, metavar="FILE",
                        help="Pelops results JSON from the matched normal sample")
    parser.add_argument("--tumor-name",  required=True, metavar="STR",
                        help="Sample identifier for the TUMOR column in the VCF")
    parser.add_argument("--normal-name", required=True, metavar="STR",
                        help="Sample identifier for the NORMAL column in the VCF")
    parser.add_argument("--output", default="-", metavar="FILE",
                        help="Output VCF path (default: stdout)")
    parser.add_argument("--srpb-threshold", type=float, default=20.0, metavar="FLOAT",
                        help="Minimum tumor SRPB to consider a rearrangement (default: 20.0)")
    parser.add_argument("--log-fold-threshold", type=float, default=3.0, metavar="FLOAT",
                        help="Minimum ln(T_SRPB/N_SRPB) for PASS on UNNAMED calls (default: 3.0, ~20x)")
    parser.add_argument("--strict-fold-change", action="store_true",
                        help="Apply fold-change filter to ALL rearrangement types including IGH-DUX4")
    parser.add_argument("--include-low-srpb", action="store_true",
                        help="Include LowSRPB records in output (default: skip them)")
    args = parser.parse_args()

    with open(args.tumor_json) as fh:
        tumor_data = json.load(fh)
    with open(args.normal_json) as fh:
        normal_data = json.load(fh)

    normal_index = load_normal_index(normal_data)

    out = open(args.output, "w") if args.output != "-" else sys.stdout

    try:
        print(
            build_header(
                normal_name=args.normal_name,
                tumor_name=args.tumor_name,
                version=tumor_data.get("version", "unknown"),
                cli_command=tumor_data.get("cli_command", ""),
                srpb_threshold=args.srpb_threshold,
                log_fold_threshold=args.log_fold_threshold,
                strict_fold_change=args.strict_fold_change,
            ),
            file=out,
        )

        for r_tumor in tumor_data.get("rearrangements", []):
            ev = r_tumor["evidence"]

            # Skip rearrangements with zero evidence in tumor
            if ev["paired_reads"] == 0 and ev["split_reads"] == 0:
                continue

            key = (r_tumor["A"]["name"], r_tumor["B"]["name"])
            r_normal = normal_index.get(key)

            records = make_bnd_pair(
                r_tumor=r_tumor,
                r_normal=r_normal,
                normal_name=args.normal_name,
                tumor_name=args.tumor_name,
                srpb_threshold=args.srpb_threshold,
                log_fold_threshold=args.log_fold_threshold,
                strict_fold_change=args.strict_fold_change,
            )

            # Optionally suppress LowSRPB records
            if not args.include_low_srpb:
                vcf_filter = records[0].split("\t")[6]
                if vcf_filter == "LowSRPB":
                    continue

            for rec in records:
                print(rec, file=out)

    finally:
        if args.output != "-":
            out.close()


if __name__ == "__main__":
    main()
