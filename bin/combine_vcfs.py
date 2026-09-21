#!/usr/bin/env python3
"""Add the germline calls (FILTER=GERMLINE, from mark_germlines.pl) to the somatic VCF of the same case.

 * samples are matched by name, so the column order of the two files does not matter
 * the germline CSQ is rewritten to the CSQ field layout of the somatic header (fields matched by
   name, fields missing in the germline VEP run are left empty, germline-only fields are dropped)
 * a variant present in both files is kept once, as the somatic record, with GERMLINE added to
   FILTER and FAIL_NVAF removed (the same rule mark_germlines.pl applies)
 * output is sorted by contig header order (natural chromosome order if not in the header) and position
"""

import argparse
import gzip
import re
import sys

HDR_KEY_RE = re.compile(r"^##(INFO|FORMAT|FILTER|ALT|contig)=<ID=([^,>]+)")
CSQ_FMT_RE = re.compile(r'^##INFO=<ID=CSQ,.*Format: ([^">]+)')
CONTIG_RE = re.compile(r"^##contig=<ID=([^,>]+)")


def read_vcf(fn):
    opener = gzip.open if fn.endswith(".gz") else open
    hdr, cols, recs = [], None, []
    with opener(fn, "rt") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            if line.startswith("##"):
                hdr.append(line)
            elif line.startswith("#CHROM"):
                cols = line.split("\t")
            else:
                recs.append(line.split("\t"))
    if cols is None:
        sys.exit(f"No #CHROM line in {fn}")
    return hdr, cols, recs


def hdr_key(line):
    m = HDR_KEY_RE.match(line)
    return f"{m.group(1)}:{m.group(2)}" if m else None


def csq_format(hdr):
    for line in hdr:
        m = CSQ_FMT_RE.match(line)
        if m:
            return m.group(1).split("|")
    return None


def build_csq_map(s_fmt, g_fmt):
    """For every somatic CSQ field, the index of the same-named germline field (-1 if missing).
    Repeated names (e.g. gnomADg_AF) are matched by occurrence."""
    g_idx = {}
    for i, name in enumerate(g_fmt):
        g_idx.setdefault(name, []).append(i)
    occ, csq_map = {}, []
    for name in s_fmt:
        n = occ.get(name, 0)
        occ[name] = n + 1
        idxs = g_idx.get(name)
        if not idxs:
            csq_map.append(-1)
        else:
            csq_map.append(idxs[n] if n < len(idxs) else idxs[-1])
    return csq_map


def remap_csq(info, csq_map):
    items = info.split(";")
    for k, item in enumerate(items):
        if not item.startswith("CSQ="):
            continue
        tx = []
        for t in item[4:].split(","):
            f = t.split("|")
            tx.append("|".join("" if i < 0 or i >= len(f) else f[i] for i in csq_map))
        items[k] = "CSQ=" + ",".join(tx)
    return ";".join(items)


def chrom_key(c):
    n = re.sub(r"^chr", "", c, flags=re.I)
    if n.isdigit():
        return (1, int(n), "")
    if n == "X":
        return (2, 0, "")
    if n == "Y":
        return (2, 1, "")
    if n.startswith("M"):
        return (2, 2, "")
    return (3, 0, n)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--somatic", required=True, help="Somatic VCF (sample order and CSQ layout of this file are kept)")
    ap.add_argument("--germline", required=True, help="Germline VCF, records flagged by mark_germlines.pl")
    args = ap.parse_args()

    s_hdr, s_cols, s_recs = read_vcf(args.somatic)
    g_hdr, g_cols, g_recs = read_vcf(args.germline)

    # ---- sample alignment ----
    s_samples = s_cols[9:]
    g_pos = {name: i for i, name in enumerate(g_cols) if i >= 9}
    for s in s_samples:
        if s not in g_pos:
            sys.exit(f"Sample {s} is in the somatic VCF but not in the germline VCF")
    if len(s_samples) != len(g_pos):
        sys.exit("Different number of samples in the somatic and germline VCF")
    g_order = list(range(9)) + [g_pos[s] for s in s_samples]

    # ---- header ----
    out_hdr = list(s_hdr)
    seen = {k for k in map(hdr_key, s_hdr) if k}
    for line in g_hdr:
        k = hdr_key(line)
        if k and k not in seen:
            out_hdr.append(line)
            seen.add(k)

    # ---- CSQ layout ----
    s_fmt, g_fmt = csq_format(s_hdr), csq_format(g_hdr)
    csq_map = build_csq_map(s_fmt, g_fmt) if s_fmt and g_fmt and s_fmt != g_fmt else []

    # ---- merge ----
    s_index = {(r[0], r[1], r[3], r[4]): i for i, r in enumerate(s_recs)}
    all_recs = list(s_recs)
    n_dup = n_new = 0
    for g in g_recs:
        rec = [g[i] for i in g_order]
        key = (rec[0], rec[1], rec[3], rec[4])
        if key in s_index:
            s = all_recs[s_index[key]]
            flt = [f for f in s[6].split(";") if f not in ("GERMLINE", "FAIL_NVAF", ".")]
            s[6] = ";".join(["GERMLINE"] + flt)
            n_dup += 1
        else:
            if csq_map:
                rec[7] = remap_csq(rec[7], csq_map)
            all_recs.append(rec)
            n_new += 1
    print(f"combine_vcfs: {n_new} germline records added, {n_dup} already present in the somatic VCF", file=sys.stderr)

    # ---- sort ----
    contigs = [m.group(1) for m in map(CONTIG_RE.match, s_hdr) if m]
    contig_rank = {c: i for i, c in enumerate(contigs)}

    def sort_key(item):
        i, r = item
        ck = (0, contig_rank[r[0]], "") if r[0] in contig_rank else chrom_key(r[0])
        return (ck, int(r[1]), i)

    out = sys.stdout
    for line in out_hdr:
        out.write(line + "\n")
    out.write("\t".join(s_cols) + "\n")
    for _, r in sorted(enumerate(all_recs), key=sort_key):
        out.write("\t".join(r) + "\n")


if __name__ == "__main__":
    main()
