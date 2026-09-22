#!/usr/bin/env python3
"""Exon regions of the genes in a mark_germlines assay json, as a bed file.

Used to restrict the germline VCF to the assay genes before VEP: the exon bed (merged, sorted, no gene
names) is overlapped with the gene intervals of the genes named in the json. Exons that only overlap
a selected gene's span are kept whole; that can include an exon of a neighbouring gene, which is harmless
because mark_germlines.pl still checks the gene symbol.

 * gene names are matched exactly against column 4 of the gene bed (split on | ; , so
   "GENE|ENSG..." style names work too)
 * genes not found in the gene bed are reported on stderr, and no match at all is an error
 * "ALL_GENES" in the json (as understood by mark_germlines.pl) keeps every exon
"""

import argparse
import bisect
import json
import re
import sys
from collections import defaultdict


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--assay", required=True, help="Assay json with a 'genes' list")
    ap.add_argument("--genes_bed", required=True, help="Gene intervals, gene name in column 4")
    ap.add_argument("--exons_bed", required=True, help="Sorted, merged exon intervals")
    args = ap.parse_args()

    with open(args.assay) as fh:
        genes = set(json.load(fh).get("genes", []))
    if not genes:
        sys.exit(f"No genes in {args.assay}")

    with open(args.exons_bed) as fh:
        exon_lines = [l.rstrip("\n").split("\t") for l in fh if l.strip() and not l.startswith(("#", "track", "browser"))]

    if "ALL_GENES" in genes:
        print("assay_gene_bed: ALL_GENES set, keeping every exon", file=sys.stderr)
        for e in exon_lines:
            print("\t".join(e))
        return

    found, spans = set(), defaultdict(list)
    with open(args.genes_bed) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            hit = genes.intersection(re.split(r"[|;,]", f[3])) if len(f) > 3 else set()
            if hit:
                found |= hit
                spans[f[0]].append((int(f[1]), int(f[2])))

    missing = sorted(genes - found)
    if not found:
        sys.exit(f"None of the {len(genes)} assay genes were found in column 4 of {args.genes_bed}")
    if missing:
        print(f"assay_gene_bed: {len(missing)} of {len(genes)} genes not found: {' '.join(missing)}", file=sys.stderr)

    # merge overlapping gene spans per contig, so an exon is tested against one sorted list
    merged = {}
    for chrom, iv in spans.items():
        iv.sort()
        out = [list(iv[0])]
        for s, e in iv[1:]:
            if s <= out[-1][1]:
                out[-1][1] = max(out[-1][1], e)
            else:
                out.append([s, e])
        merged[chrom] = ([s for s, _ in out], [e for _, e in out])

    n = 0
    for e in exon_lines:
        m = merged.get(e[0])
        if not m:
            continue
        starts, ends = m
        s, t = int(e[1]), int(e[2])
        i = bisect.bisect_left(starts, t) - 1   # last span starting before the exon ends
        if i >= 0 and ends[i] > s:
            print("\t".join(e))
            n += 1
    print(f"assay_gene_bed: {len(found)} genes, {n} exon intervals", file=sys.stderr)


if __name__ == "__main__":
    main()
