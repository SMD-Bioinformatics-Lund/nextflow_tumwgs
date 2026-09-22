# Known Issues

## Solid assay gene list uses stale HGNC symbols, 9 genes silently drop out of germline reporting

**Found:** 2026-09-22, while restricting germline `INTERSECT_CODING` to the assay gene panel
(`bin/assay_gene_bed.py`, `PANEL_EXON_BED` in `modules/local/filters/main.nf`).

**What happens:** `resources/solid_mg.json` lists 210 gene symbols. When matched against
`gencode.v33.annotation.genes.proteincoding.bed` (column 4, the source `params.gencode_genes`
points at), only 200 are found. The other 10 are:

```
C19MC GPR126 H2AFX H3F3A H3F3B HIST1H3B HIST1H3C HRPT2 MRE11A TCEB1
```

- `C19MC` is a non-coding miRNA cluster — correctly absent from a protein-coding gene bed, not a
  bug.
- The remaining 9 look like older HGNC symbols that Gencode v33 has since renamed:

  | solid_mg.json | current HGNC symbol (Gencode v33) |
  |---|---|
  | GPR126 | ADGRG6 |
  | H2AFX | H2AX |
  | H3F3A | H3-3A |
  | H3F3B | H3-3B |
  | HIST1H3B | H3C2 |
  | HIST1H3C | H3C3 |
  | HRPT2 | CDC73 |
  | MRE11A | MRE11 |
  | TCEB1 | ELOC |

**Impact:** `assay_gene_bed.py` warns on stderr and continues (by design — one bad gene name
shouldn't fail the run), but the practical effect is that these 9 genes get **no exon regions in
the germline panel bed**, so any germline variant in them is invisible to `MARK_GERMLINES` for the
`solid` profile. This predates the `PANEL_EXON_BED` change — `mark_germlines.pl` matches on the
same symbols against VEP's `SYMBOL` field, so these genes were already silently excluded from
germline calling; restricting the exon bed just surfaced it via the stderr mismatch report.

**Fix options (not yet applied):**
1. Update `resources/solid_mg.json` to the current HGNC symbols above.
2. Have `assay_gene_bed.py` (and/or `mark_germlines.pl`) accept a small alias table for known
   renames, so older assay JSONs keep working without edits.

Option 1 is simplest and should be done with clinical sign-off on the gene list, same as the
open hema gene-list item.
