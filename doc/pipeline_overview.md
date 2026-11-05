# Somatic Whole Genome Sequencing (SWGS) Pipeline — Overview

**Organisation:** SMD Bioinformatics, Lund
**Pipeline version:** 3.0.0
**Reference genome:** GRCh38

---

## What is this pipeline?

The SWGS pipeline analyses tumour DNA from whole-genome sequencing (WGS) to detect cancer-relevant genetic variants. Given raw sequencing data from a tumour biopsy and a matched normal (blood) sample, the pipeline:

1. Aligns reads to the human reference genome
2. Detects somatic mutations — changes present in the tumour but not in the normal
3. Annotates variants with clinical and population databases
4. Filters results to a clinically relevant gene panel
5. Loads findings into a clinical review database (Coyote)

The pipeline is designed for clinical-grade analysis at a haematology and oncology centre and runs on high-performance computing (HPC) clusters.

---

## Who is this for?

This pipeline is used by clinical bioinformaticians and laboratory geneticists at SMD Bioinformatics, Lund. Results feed into clinical decision-making for patients with haematologic malignancies and solid tumours.

---

## What does the pipeline detect?

### Single Nucleotide Variants (SNVs) and small insertions/deletions (indels)

Point mutations and small insertions or deletions within protein-coding genes. These are the most common class of somatic driver mutations in cancer.

Four independent callers are run and their results are combined, increasing sensitivity while allowing cross-validation:

| Caller | Approach |
|--------|----------|
| **Sentieon TNscope** | Machine learning-based somatic variant detection |
| **DeepSomatic** | Deep neural network (derived from Google DeepVariant) |
| **FreeBayes** | Bayesian haplotype-based variant detection |
| **VarDict** | Allele frequency-based somatic variant detection |
| **Pindel** | Specialised caller for small-to-medium insertions and deletions |

Additionally, **germline variants** are called from the normal sample using Sentieon DNAscope.

### Copy Number Variants (CNVs)

Gains and losses of genomic segments — for example, amplification of an oncogene or deletion of a tumour suppressor gene. CNVs are detected using GATK4's read-depth and allele-frequency-based method.

### Loss of heterozygosity and chromosome count (tumour–normal pairs only)

Every chromosome normally carries one copy from each parent. **Loss of heterozygosity (LOH)** means the tumour has lost one parent's copy in a region, either with a loss of DNA or with the remaining copy duplicated (copy-neutral LOH, which a read-depth CNV caller cannot see). LOH can uncover a mutated tumour suppressor gene and is part of several diagnostic classifications.

For each tumour–normal pair, the pipeline uses **ASCAT** to count the copies of each parental chromosome copy along the genome. From that it reports:

| Result | What it means |
|---|---|
| Tumour purity and ploidy | The share of tumour cells in the sample, and the average number of chromosome copies per tumour cell. When the tumour genome has too few changes, purity cannot be measured and is reported as such |
| LOH regions | Region, type (with a loss, copy-neutral, or with a gain) and the share of tumour cells that carry it |
| Chromosome count | The number of chromosomes in the tumour, counted the way a karyotype counts them: a chromosome is gained or lost only when the whole chromosome is, while changes of one arm (for example del(7q)) are listed separately. A range is given when a chromosome is on the edge |
| Gains of part of an arm | Gains covering part of a chromosome arm (at least 5 Mb), with the share of tumour cells that carry them. They are listed, not counted |
| High hyperdiploidy (HeH) | For acute lymphoblastic leukaemia: 51–67 chromosomes, with a check for doubled hypodiploidy, which can reach the same count. Not given for solid tumours |
| Sample checks | Whether the sex in the data matches the sample sheet, and the share of tumour cells found in the normal sample |
| Changes in part of the tumour | Gains and LOH present in only some tumour cells (subclonal), with an estimate of how many |

Limitations:
- Changes present in fewer than about 12–15% of tumour cells are not detected.
- Loss of the Y chromosome is not assessed.
- The thresholds have been compared against another tool (PURPLE) on a small set of samples but not yet validated against karyotype or FISH.
- Results are written as files; they are not yet shown in Coyote.

### Structural Variants (SVs) and Gene Fusions

Large-scale rearrangements of the genome, including translocations, large deletions, inversions, and duplications. When a rearrangement joins two genes together, it can create an oncogenic **fusion gene** (e.g., BCR-ABL1 in CML). SVs are detected using **Manta** and annotated using **SnpEff**.

---

## Pipeline inputs

The pipeline requires:

- **Paired-end FASTQ files** — raw sequencing reads from Illumina whole-genome sequencing (typically 30–80× coverage)
- **A matched normal sample** — blood or buccal DNA from the same patient, used to distinguish somatic (tumour-specific) from germline (inherited) changes
- **A sample description CSV** — a file listing sample IDs, file paths, diagnosis, and optional metadata (tumour purity, FFPE status, etc.)

### Input CSV format

```
id,group,diagnosis,type,read1,read2,...
TUMOR_01,GRP01,AML,T,/data/TUMOR_01_R1.fastq.gz,/data/TUMOR_01_R2.fastq.gz,...
NORMAL_01,GRP01,AML,N,/data/NORMAL_01_R1.fastq.gz,/data/NORMAL_01_R2.fastq.gz,...
```

Tumour and normal samples are paired by the `group` field. The `type` field must be `T` (tumour) or `N` (normal).

---

## Disease profiles

The pipeline can be run in two modes, each using a disease-appropriate gene panel:

| Profile | Disease context | Gene panel scope |
|---------|----------------|-----------------|
| **hema** | Haematologic malignancies (leukaemia, lymphoma, myeloma) | Haematology-specific SNV, CNV, and fusion gene panels |
| **solid** | Solid tumours (lung, colorectal, breast, etc.) | BTB solid tumour SNV, CNV, and fusion gene panels |

Variant calls outside the relevant gene panel are filtered out before clinical reporting.

---

## Analysis stages

```
Raw FASTQ
    │
    ▼
Quality control & optional trimming
    │
    ▼
Alignment to GRCh38 (Sentieon BWA-MEM)
    │
    ├── Duplicate marking
    ├── Base quality score recalibration (BQSR)
    └── Alignment quality metrics
         │
         ▼
    ┌────┴───────────────────────────┐
    │                                │
    │    Sample identity check       │
    │    (ID-SNP verification)       │
    │                                │
    └────────────────────────────────┘
         │
         ▼
    ┌──────────────────────────────────────────────────────────┐
    │ Parallel variant calling                                  │
    │                                                           │
    │  SNVs/indels         CNVs              SVs & fusions      │
    │  (5 callers)         (GATK4)           (Manta + SnpEff)   │
    │                      LOH, ploidy,                         │
    │                      chromosome count                     │
    │                      (ASCAT, pairs)                       │
    └──────────────────────────────────────────────────────────┘
         │
         ▼
Variant annotation
  · Consequence prediction (Ensembl VEP v113)
  · Population frequencies (gnomAD v4.0)
  · Pathogenicity scoring (CADD v1.7)
  · Cancer mutation database (COSMIC v92)
  · Germline variant identification
         │
         ▼
Gene panel filtering
(Disease-specific SNV, CNV, fusion panels)
         │
         ▼
Panel-of-normals filtering
(Remove recurrent sequencing artefacts)
         │
         ▼
Clinical database import (Coyote)
```

---

## Pipeline outputs

All results are written to the output directory specified at runtime.

### Variant calls

| Output | Format | Description |
|--------|--------|-------------|
| Final SNV/indel calls | VCF | Annotated, filtered somatic variants in panel genes |
| Germline variants | VCF | Inherited variants from normal sample (DNAscope) |
| CNV segments | BED / VCF | Copy number gains and losses with gene annotations |
| LOH and chromosome count | TXT / TSV / JSON | Purity, ploidy, LOH regions, gains, chromosome count, HeH verdict, sample checks (tumour–normal pairs) |
| Fusion gene calls | VCF / TSV | Structural variants predicted to create gene fusions |

### Quality metrics

| Output | Description |
|--------|-------------|
| Alignment QC | Coverage, mapping rate, insert size distribution |
| Duplication rate | Proportion of PCR duplicates in the library |
| Identity verification | Confirms tumour and normal are from the same patient |

### Visualisations

| Output | Description |
|--------|-------------|
| CNV plot | Genome-wide copy number profile |
| Coverage plot | Read depth across the genome |

### Database export

Results are automatically imported into **Coyote**, the clinical variant review database used for patient reporting. A structured YAML file is also generated for archiving.

### Run metadata

| Output | Description |
|--------|-------------|
| `software_versions.yml` | Exact version of every tool used in the run |
| Nextflow report (HTML) | Process timing, resource usage, success/failure summary |
| Nextflow trace (TXT) | Per-process execution log |

---

## Annotation databases

Variants are annotated against the following databases:

| Database | Version | Purpose |
|----------|---------|---------|
| Ensembl VEP | v113.0 | Gene consequence, HGVS nomenclature |
| gnomAD | v4.0 | Population allele frequencies (germline filtering) |
| CADD | v1.7 | Computational pathogenicity score |
| COSMIC | v92 | Curated somatic cancer mutations |
| dbSNP | v146 | Known variant identifiers |
| ClinVar | — | Clinical significance classifications |

---

## Infrastructure & reproducibility

- **Compute:** Runs on SLURM HPC clusters (Grace/Hopper at Lund)
- **Containers:** All tools run inside a Singularity container (`tumwgs_container.sif`), ensuring full reproducibility across runs and cluster environments
- **Workflow manager:** Nextflow DSL2 — provides automatic parallelisation, failure recovery, and run logging
- **Resume capability:** Failed or interrupted runs can be resumed from the last successful step without re-running the full pipeline

---

## Version history

See [CHANGELOG.md](../CHANGELOG.md) for a full list of changes between pipeline versions.

---

## Contact

For questions about running the pipeline or interpreting results, contact the SMD Bioinformatics team at Klinisk genetik patologi och molekylär diagnostik at Region Skåne.
For bug reports and feature requests, open an issue on the [GitHub repository](https://github.com/SMD-Bioinformatics-Lund/nextflow_tumwgs).
