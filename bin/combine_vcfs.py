#!/usr/bin/env python3
"""Add the germline calls (FILTER=GERMLINE, from mark_germlines.pl) to the somatic VCF of the same case.

 * samples are matched by name, so the column order of the two files does not matter
 * the germline CSQ is rewritten to the CSQ field layout of the somatic header (fields matched by
   name, fields missing in the germline VEP run are left empty, germline-only fields are dropped)
 * a variant present in both files is kept once, as the somatic record, with GERMLINE added to
   FILTER and FAIL_NVAF removed (the same rule mark_germlines.pl applies)
 * output is sorted by contig header order (natural chromosome order if not in the header) and position

The germline VCF only holds the flagged variants, so it is kept in memory; the somatic VCF is streamed
(one scan to check the sort order and find the shared variants, one to merge), so memory does not grow
with the number of somatic variants. A somatic VCF that is not sorted is sorted in memory instead.
"""

import argparse
import gzip
import re
import sys

HDR_KEY_RE = re.compile(r"^##(INFO|FORMAT|FILTER|ALT|contig)=<ID=([^,>]+)")
CSQ_FMT_RE = re.compile(r'^##INFO=<ID=CSQ,.*Format: ([^">]+)')
CONTIG_RE = re.compile(r"^##contig=<ID=([^,>]+)")


def open_text(fn):
    return gzip.open(fn, "rt") if fn.endswith(".gz") else open(fn, "rt")


def read_header(fh, fn):
    """Read the ## lines and the #CHROM line, leaving fh at the first record."""
    hdr, cols = [], None
    for line in fh:
        line = line.rstrip("\n")
        if not line.strip():
            continue
        if line.startswith("##"):
            hdr.append(line)
        elif line.startswith("#CHROM"):
            cols = line.split("\t")
            break
        else:
            break
    if cols is None:
        sys.exit(f"No #CHROM line in {fn}")
    return hdr, cols


def read_vcf(fn):
    with open_text(fn) as fh:
        hdr, cols = read_header(fh, fn)
        recs = [line.rstrip("\n").split("\t") for line in fh if line.strip()]
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

    g_hdr, g_cols, g_recs = read_vcf(args.germline)
    with open_text(args.somatic) as fh:
        s_hdr, s_cols = read_header(fh, args.somatic)

        # ---- sample alignment ----
        s_samples = s_cols[9:]
        g_pos = {name: i for i, name in enumerate(g_cols) if i >= 9}
        for s in s_samples:
            if s not in g_pos:
                sys.exit(f"Sample {s} is in the somatic VCF but not in the germline VCF")
        if len(s_samples) != len(g_pos):
            sys.exit("Different number of samples in the somatic and germline VCF")
        g_order = list(range(9)) + [g_pos[s] for s in s_samples]
        g_recs = [[g[i] for i in g_order] for g in g_recs]

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

        # ---- sort key: contig header order, else natural chromosome order; then position ----
        contig_rank = {m.group(1): i for i, m in enumerate(filter(None, map(CONTIG_RE.match, s_hdr)))}
        ck_cache = {}

        def sort_key(chrom, pos):
            ck = ck_cache.get(chrom)
            if ck is None:
                ck = ck_cache[chrom] = (0, contig_rank[chrom], "") if chrom in contig_rank else chrom_key(chrom)
            return (ck, int(pos))

        # ---- pass 1: is the somatic VCF sorted, and which germline variants does it already hold? ----
        g_keys = {(g[0], g[1], g[3], g[4]) for g in g_recs}
        g_sites = {(k[0], k[1]) for k in g_keys}
        present, is_sorted, prev = set(), True, None
        for line in fh:
            if not line.strip():
                continue
            f = line.split("\t", 5)
            sk = sort_key(f[0], f[1])
            if prev is not None and sk < prev:
                is_sorted = False
            prev = sk
            if (f[0], f[1]) in g_sites:
                key = (f[0], f[1], f[3], f[4])
                if key in g_keys:
                    present.add(key)

    n_dup = sum(1 for g in g_recs if (g[0], g[1], g[3], g[4]) in present)
    new_germ = [g for g in g_recs if (g[0], g[1], g[3], g[4]) not in present]
    if csq_map:
        for g in new_germ:
            g[7] = remap_csq(g[7], csq_map)
    new_germ = sorted(((sort_key(g[0], g[1]), i, g) for i, g in enumerate(new_germ)), key=lambda t: t[:2])
    print(f"combine_vcfs: {len(new_germ)} germline records added, {n_dup} already present in the somatic VCF", file=sys.stderr)

    def add_germline_flag(line):
        r = line.rstrip("\n").split("\t")
        if (r[0], r[1], r[3], r[4]) not in present:
            return line
        flt = [f for f in r[6].split(";") if f not in ("GERMLINE", "FAIL_NVAF", ".")]
        r[6] = ";".join(["GERMLINE"] + flt)
        return "\t".join(r) + "\n"

    # ---- pass 2: merge the sorted germline records into the somatic stream ----
    out = sys.stdout
    for line in out_hdr:
        out.write(line + "\n")
    out.write("\t".join(s_cols) + "\n")

    with open_text(args.somatic) as fh:
        read_header(fh, args.somatic)
        recs = (line if line.endswith("\n") else line + "\n" for line in fh if line.strip())
        if is_sorted:
            stream = ((sort_key(*l.split("\t", 2)[:2]), l) for l in recs)
        else:
            print("combine_vcfs: somatic VCF is not sorted, sorting in memory", file=sys.stderr)
            keyed = [(sort_key(*l.split("\t", 2)[:2]), i, l) for i, l in enumerate(recs)]
            keyed.sort(key=lambda t: t[:2])
            stream = ((k, l) for k, _, l in keyed)

        gi = 0
        for sk, line in stream:
            # germline-only records go before a somatic record only when strictly earlier, so
            # somatic records stay ahead of germline records at the same position
            while gi < len(new_germ) and new_germ[gi][0] < sk:
                out.write("\t".join(new_germ[gi][2]) + "\n")
                gi += 1
            out.write(add_germline_flag(line) if present else line)
        for _, _, g in new_germ[gi:]:
            out.write("\t".join(g) + "\n")


if __name__ == "__main__":
    main()
