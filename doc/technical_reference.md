# SWGS Pipeline — Technical Reference

**Pipeline:** Somatic Whole Genome Sequencing (SWGS)
**Version:** 3.0.0
**Framework:** Nextflow DSL2
**Genome build:** GRCh38 (GCA_000001405.15, no-alt)

---

## Table of Contents

1. [Architecture Overview](#1-architecture-overview)
2. [Directory Layout](#2-directory-layout)
3. [Input Specification](#3-input-specification)
4. [Running the Pipeline](#4-running-the-pipeline)
5. [Configuration System](#5-configuration-system)
6. [Subworkflows](#6-subworkflows)
   - 6.1 [Input Validation & Metadata (CHECK_INPUT)](#61-check_input)
   - 6.2 [Downsampling & Trimming (SAMPLE)](#62-sample)
   - 6.3 [Alignment (ALIGN_SENTIEON)](#63-align_sentieon)
   - 6.4 [Quality Control (BAM_QC)](#64-bam_qc)
   - 6.5 [SNV Calling (SNV_CALLING)](#65-snv_calling)
   - 6.6 [SNV Annotation (SNV_ANNOTATE)](#66-snv_annotate)
   - 6.7 [CNV Calling (CNV_CALLING_WGS)](#67-cnv_calling_wgs)
   - 6.8 [SV Calling (SV_CALLING)](#68-sv_calling)
   - 6.9 [Visualisation (VISUALIZE)](#69-visualize)
   - 6.10 [Database Import (ADD_TO_DB)](#610-add_to_db)
7. [Profiles](#7-profiles)
8. [Reference Data & Annotation Databases](#8-reference-data--annotation-databases)
9. [Software Versions](#9-software-versions)
10. [Output Structure](#10-output-structure)
11. [Adding New Components](#11-adding-new-components)
12. [Parameter Reference](#12-parameter-reference)

---

## 1. Architecture Overview

The pipeline is structured in three layers:

```
main.nf
  └── workflows/common.nf          (SWGP_COMMON — full DAG)
        ├── subworkflows/local/     (logical stage groupings)
        └── modules/local/          (individual process definitions)
```

**Execution flow:**

```
CSV ──► CHECK_INPUT ──► SAMPLE ──► ALIGN_SENTIEON
                                        │
                                   ┌────┴────┐
                                BAM_QC    (parallel variant calling)
                                        │
                          ┌─────────────┼─────────────┐
                     SNV_CALLING  CNV_CALLING_WGS  SV_CALLING
                          │
                     SNV_ANNOTATE
                          │
                       VISUALIZE
                          │
                       ADD_TO_DB
                          │
              CUSTOM_DUMPSOFTWAREVERSIONS
```

All variant-calling stages run in parallel after alignment is complete.

---

## 2. Directory Layout

```
nextflow_tumwgs/
├── main.nf                         # Entry point; delegates to workflows/common.nf
├── nextflow.config                 # Primary configuration (profiles, params, containers)
├── README.MD                       # Quick-start guide
├── CLAUDE.md                       # Developer guide and documentation plan
├── CHANGELOG.md                    # Version history
├── Makefile                        # Build targets (docs, linting, etc.)
├── requirements.txt                # Python dependencies for helper scripts
│
├── workflows/
│   └── common.nf                   # SWGP_COMMON: full pipeline DAG
│
├── subworkflows/local/
│   ├── create_meta.nf              # CHECK_INPUT
│   ├── sample.nf                   # SAMPLE (downsampling/trimming)
│   ├── align_sentieon.nf           # ALIGN_SENTIEON
│   ├── bam_qc.nf                   # BAM_QC
│   ├── snv_calling.nf              # SNV_CALLING
│   ├── snv_annotate.nf             # SNV_ANNOTATE
│   ├── cnv_calling_wgs.nf          # CNV_CALLING_WGS
│   ├── sv_calling.nf               # SV_CALLING
│   ├── visualize.nf                # VISUALIZE
│   ├── add_to_db.nf                # ADD_TO_DB
│   ├── check_idsnp.nf              # Identity SNP verification
│   ├── cnv_annotate.nf             # CNV annotation helper
│   └── cnv_calling.nf              # Panel CNV (legacy/panel-mode)
│
├── modules/local/
│   ├── GATK/                       # GATK4 CNV processes
│   ├── ascat/                      # ASCAT allele-specific CN / LOH / HeH
│   ├── deepSomatic/                # DeepSomatic variant caller
│   ├── fastp/                      # FASTQ trimming
│   ├── filters/                    # VCF filtering scripts
│   ├── freebayes/                  # Freebayes SNV caller
│   ├── gens/                       # GENS visualisation
│   ├── idSnp/                      # ID-SNP genotyping
│   ├── manta/                      # Manta SV caller
│   ├── pindel/                     # Pindel indel caller
│   ├── qc/                         # QC metric collection
│   ├── sentieon/                   # Sentieon alignment + variant calling
│   ├── seqtk/                      # Downsampling
│   ├── snpeff/                     # SnpEff SV annotation
│   ├── vardict/                    # VarDict SNV caller
│   ├── concatenate_vcfs/           # VCF merging
│   ├── coyote/                     # DB import
│   ├── custom/                     # Miscellaneous helpers
│   └── check_input/                # CSV validation
│
├── configs/
│   ├── modules/                    # Per-subworkflow resource configs
│   │   ├── base.config             # Default labels (process_high, etc.)
│   │   ├── align_sentieon.config
│   │   ├── snv_calling.config
│   │   ├── snv_annotate.config
│   │   ├── cnv_calling.config
│   │   ├── sv_calling.config
│   │   └── ...
│   └── nextflow.*.config           # Cluster-specific configs
│
├── bin/                            # 46 helper scripts (Python, Perl, Bash, R)
├── resources/                      # Shards, gene lists, interval BED files
└── doc/                            # This documentation
```

---

## 3. Input Specification

### 3.1 CSV Format

The pipeline takes a CSV file via `--csv`. Each row represents one sample (fastq pair). Tumor and matched normal are linked by a shared `group` value.

```csv
id,group,diagnosis,type,clarity_sample_id,read1,read2,ffpe,purity,sequencing_run,n_reads,clarity_pool_id
TUMOR1,GRP01,AML,T,CLR001,/data/TUMOR1_R1.fastq.gz,/data/TUMOR1_R2.fastq.gz,false,0.75,RUN001,80000000,POOL01
NORMAL1,GRP01,AML,N,CLR002,/data/NORMAL1_R1.fastq.gz,/data/NORMAL1_R2.fastq.gz,false,1.0,RUN001,80000000,POOL01
```

**Required columns:**

| Column | Type | Description |
|--------|------|-------------|
| `id` | string | Unique sample identifier |
| `group` | string | Links tumor and matched normal together |
| `diagnosis` | string | Clinical diagnosis code |
| `type` | `T` or `N` | Tumor or normal |
| `clarity_sample_id` | string | LIMS sample ID |
| `read1` | path | Absolute path to R1 FASTQ (gzip) |
| `read2` | path | Absolute path to R2 FASTQ (gzip) |

**Optional columns:**

| Column | Default | Description |
|--------|---------|-------------|
| `ffpe` | `false` | FFPE sample — enables FFPE-aware variant filters |
| `purity` | `""` | Tumor purity estimate (0–1) |
| `sequencing_run` | `""` | Run ID for traceability |
| `n_reads` | `""` | Target read count for downsampling (leave empty to disable) |
| `clarity_pool_id` | `""` | Library pool ID |

### 3.2 Channel Schema

After `CHECK_INPUT`, two channels are created:

- `fastq_ch`: `[group, meta, fastq_R1, fastq_R2]`
- `meta_ch`: `[group, meta]`

The `meta` map carries all CSV fields plus derived fields (e.g., `meta.sex`, `meta.sub`, `meta.tumor_only`).

### 3.3 Tumor-only Mode

A `group` may contain either two rows (tumor + matched normal, `type` = `T`/`N`) or a single tumor row. `check_samplesheet.py` validates this at input time — any group with 0, 1-normal-only, 2-tumors, or 3+ rows is rejected before the pipeline runs. `create_meta.nf` derives `meta.tumor_only` (`true` when the group has a single row) and stamps it on every sample's meta map, so downstream logic can branch on an explicit flag rather than re-deriving pairing from list sizes.

Sample count per group (1 vs 2) is still the mechanism most processes branch on internally (`meta.id.size()`), since `meta` becomes a list once grouped by `group`. Tumor-only support differs by stage:

| Stage | Tumor-only behavior |
|-------|----------------------|
| SNV calling (Freebayes, VarDict, TNscope, DeepSomatic, Pindel) | Each caller runs a native single-sample mode (DeepSomatic uses its dedicated `WGS_TUMOR_ONLY` model). No matched-normal comparison. |
| DNAscope | Unaffected — always germline, per-sample regardless of pairing. |
| PON_FILTER | Unaffected — already only references the tumor sample. |
| Germline/somatic filtering | No matched-normal genotype subtraction is available. Somatic specificity relies entirely on the Panel of Normals (`PON_freebayes`/`PON_vardict`) plus population allele frequency filtering (`POST_ANNOTATION_FILTERS`, gnomAD `MAX_AF` vs `params.filter_freq`). Expect a higher germline/artifact carry-through than paired calling. |
| GATK CNV | `GATKCOV_CALL_GERMLINE` (segment calling on the normal) is skipped for tumor-only groups — there is no normal to call. `GATKCOV_CALL`'s `--normal-allelic-counts` refinement is likewise omitted. Sex-specific PON (`GATK_PON_FEMALE`/`GATK_PON_MALE`) substitutes for the normal at the read-count level regardless of pairing. |
| ASCAT (LOH/HeH) | Skipped — ASCAT needs the matched normal's heterozygous SNPs. Paired groups are also skipped (with a `[ASCAT]` warning) when the normal has no `F`/`M` sex. |
| Manta | Runs `configManta.py --tumorBam` only (no `--normalBam`); emits `tumorSV.vcf.gz` as the sole output. The `manta_vcf_normal` channel is declared `optional: true` and is simply absent for these groups. |
| QC (`SOMALIER_QC`) | Runs single-sample somalier extraction/relate; the pairwise contamination/pairedness comparison is skipped since there's no normal to compare against. `VERIFYBAMID2` (population-based contamination) and the ID-SNP checks (`SNP_CHECK`, `PAIRGEN_CDM`) already run per-sample and report `is_paired_sample: false`. |

---

## 4. Running the Pipeline

### 4.1 Manual Execution

```bash
module load singularity Java nextflow/21.10.6

nextflow run main.nf \
    -entry SWGP \
    -c nextflow.config \
    --csv /path/to/samples.csv \
    -profile hema \
    -with-report /path/to/reports/sample.report.html \
    -with-trace  /path/to/reports/sample.trace.txt \
    -with-timeline /path/to/reports/sample.timeline.html \
    -work-dir /path/to/work
```

Replace `-profile hema` with `-profile solid` for solid tumour samples.

### 4.2 Automated Production Run

In production, execution is triggered by `start_nextflow_analysis.pl` (from the `bnf-infrastructure` repository). It monitors for input CSV files produced by the `Bjorn` LIMS system and triggers the pipeline automatically.

Pipeline identity in `pipeline_files.config`:

```ini
[tumwgs-hema]
pipeline = /production/nextflow_tumwgs/main.nf -entry SWGP --profile hema
container = /production/nextflow_tumwgs/container/tumwgs_container.sif
singularity_version = 3.8.0
nextflow_version = 21.04.2
executor = slurm
cluster = grace
queue = normal
```

### 4.3 Development / Test Run

```bash
nextflow run main.nf -entry SWGP -c nextflow.config --csv sample.csv \
    -profile test --dev
```

The `test` profile uses reduced gene panels and lower resource allocations. Disable database import during testing with `--coyote false`.

### 4.4 Resume

```bash
nextflow run main.nf ... -resume
```

Nextflow caches all completed processes. `-resume` re-uses cached results and only re-runs processes whose inputs have changed.

---

## 5. Configuration System

### 5.1 Primary Config (`nextflow.config`)

Contains:
- Container path (`process.container`)
- Sentieon licence server (`env.SENTIEON_LICENSE`)
- SLURM executor settings
- Global `params` defaults
- Profile blocks for `hema`, `solid`, `test`, `hopper`, `trannel`

### 5.2 Module Resource Configs (`configs/modules/*.config`)

Each subworkflow has a dedicated resource config included from `nextflow.config`. Labels map to CPU/memory/time:

| Label | CPUs | Memory | Wall time |
|-------|------|--------|-----------|
| `process_high` | 50 | varies | 48 h |
| `process_medium` | 16 | varies | 48 h |
| `process_low` | 8 | varies | 48 h |
| `process_single` | 1 | varies | 48 h |

Override per-process in the module's config block.

### 5.3 Profile Resolution

Profiles stack. A typical production run uses two profiles:

```
-profile hopper,hema
```

- `hopper` sets the executor and container
- `hema` sets the gene panels and coyote group

---

## 6. Subworkflows

### 6.1 CHECK_INPUT

**File:** `subworkflows/local/create_meta.nf`

Validates the input CSV and emits sample channels.

**Key functions:**
- `CSV_CHECK` — header validation, required field checks
- `create_fastq_channel()` — builds `[group, meta, R1, R2]` tuples
- `create_samples_channel()` — builds `[group, meta]` tuples

**Outputs:**
- `fastq` — per-sample FASTQ channel
- `meta` — per-sample metadata channel

---

### 6.2 SAMPLE

**File:** `subworkflows/local/sample.nf`

Optional downsampling and adapter trimming before alignment.

**Processes:**

| Process | Tool | Trigger |
|---------|------|---------|
| `SEQTK` | seqtk sample | `meta.sub != false` (i.e., `n_reads` set in CSV) |
| `FASTP` | fastp | `params.trimfq == true` |

**Output:** Processed FASTQ (or passthrough if neither step enabled)

---

### 6.3 ALIGN_SENTIEON

**File:** `subworkflows/local/align_sentieon.nf`

Sharded BWA-MEM alignment via Sentieon, followed by deduplication and BQSR.

**Process chain:**

```
BWA_ALIGN_SHARD ×8 ──► BWA_MERGE_SHARDS ──► BAM_CRAM
                                                  │
                                             MARKDUP
                                                  │
                                        REALIGN_INDEL_BQSR
                                                  │
                                           CRAM_TO_BAM
```

**Key parameters:**

| Param | Default | Description |
|-------|---------|-------------|
| `params.bwa_shards` | `8` | Number of parallel alignment shards |
| `params.K_size` | `100000000` | Reads per processing block |
| `params.known` | — | Known indels VCF (Mills + 1000G) |
| `params.dbSnp` | — | dbSNP VCF for BQSR |

**Outputs:**
- `bam_bqsr` — Final analysis-ready BAM (BQSR applied)
- `cram_bqsr` — CRAM equivalent
- `cram_dedup` — Deduplicated CRAM (pre-BQSR, for archiving)
- `dedup_metrics` — Sentieon duplicate metrics

---

### 6.4 BAM_QC

**File:** `subworkflows/local/bam_qc.nf`

Alignment quality metrics and sample identity verification.

**Processes:**

| Process | Tool | Output |
|---------|------|--------|
| `SENTIEON_QC` | Sentieon QualityMap | Per-base and aggregate alignment stats |
| `COLLECT_QC` | Custom Python | Summary QC table |
| `QC_TO_CDM` | Custom | CDM-formatted QC |
| `ALLELE_CALL` | Custom | Genotypes at ID-SNP loci |
| `SNP_CHECK` | Custom | Tumor/normal identity comparison |
| `PAIRGEN_CDM` | Custom | Pairwise comparison export |

The ID-SNP check genotypes samples at ~50 pre-defined SNP positions to confirm that the tumour and normal are from the same individual. For tumor-only groups, `SOMALIER_QC` skips the pairwise comparison entirely — see [§3.3 Tumor-only Mode](#33-tumor-only-mode).

---

### 6.5 SNV_CALLING

**File:** `subworkflows/local/snv_calling.nf`

Somatic SNV and small indel calling using four callers. Calls are made per genomic region (BED intervals) and then merged.

**Callers:**

| Process | Tool | Mode | Notes |
|---------|------|------|-------|
| `FREEBAYES` | FreeBayes | Haplotype-based | Bayesian variant calling |
| `VARDICT` | VarDict Java | Somatic | Min VAF: `params.vardict_var_freq_cutoff_p` |
| `TNSCOPE_ML` | Sentieon TNscope | Tumor/normal ML | Min VAF: `params.tnscope_var_freq_cutoff_p` |
| `DEEPSOMATIC` | DeepSomatic v1.9.0 | Deep learning | GPU-accelerated when available |
| `PINDEL_CALLING` | Pindel | Indels | Split/read-pair indel calling |
| `DNASCOPE` | Sentieon DNAscope | Germline | Germline variants from normal |

**Aggregation:**
1. Per-caller, per-region VCFs → `CONCATENATE_VCFS` (per caller)
2. All callers → `AGGREGATE_VCFS` (combined VCF with caller support tags)

**VAF thresholds (default):**

| Caller | Paired VAF | Tumor-only VAF |
|--------|-----------|----------------|
| Freebayes | 0.03 | 0.05 |
| VarDict | 0.03 | 0.05 |
| TNscope | 0.01 | 0.05 |

**Toggle any caller** in `nextflow.config`:
```groovy
params.freebayes   = true
params.vardict     = true
params.tnscope     = true
params.deepsomatic = true
params.pindel      = true
params.dnascope    = true
```

---

### 6.6 SNV_ANNOTATE

**File:** `subworkflows/local/snv_annotate.nf`

Annotation and filtering of aggregated SNV calls.

**Process chain:**

```
agg_vcf ──► PON_FILTER ──► ANNOTATE_VEP ──► FILTER_PANEL ──► FIX_VEP ──► POST_ANNOTATION_FILTERS
```

**Processes:**

| Process | Purpose |
|---------|---------|
| `PON_FILTER` | Remove recurrent artefacts using Panel of Normals (206 normals) |
| `ANNOTATE_VEP` | Ensembl VEP v113 — consequence, CADD, gnomAD, COSMIC, HGVS |
| `FILTER_PANEL` | Keep only variants in `params.PANEL_SNV` gene list |
| `FIX_VEP` | Reformat VEP CSQ field for downstream compatibility |
| `POST_ANNOTATION_FILTERS` | Filter by population frequency, FILTER field, pathogenicity override |

**Key filter parameters:**

| Param | Default | Description |
|-------|---------|-------------|
| `params.filter_freq` | `0.05` | Maximum gnomAD allele frequency |
| `params.filter_field_filter` | `"FAIL_PON*,FAIL_NVAF,FAIL_LONGDEL"` | FILTER field values to exclude |
| `params.override_filter_terms` | `"CLIN_SIG=likely_pathogenic,pathogenic"` | Override filter for known pathogenic |

**VEP annotation sources:**

| Database | Version | Field |
|----------|---------|-------|
| gnomAD | v4.0 | AF_gnomAD |
| CADD | v1.7 | CADD_PHRED |
| COSMIC | v92 | COSMIC_ID |
| dbSNP | v146 | Existing_variation |
| ClinVar | — | CLIN_SIG |

---

### 6.7 CNV_CALLING_WGS

**File:** `subworkflows/local/cnv_calling_wgs.nf`

GATK4-based somatic copy-number variant calling for WGS data.

**Process chain:**

```
bam ──► GATKCOV_COUNT ──┐
                         ├──► GATKCOV_CALL ──► OVERLAP_GENES ──► FILTER_CNVS_PANEL
bam ──► GATKCOV_BAF ───┘
                         └──► GATKCOV_CALL_GERMLINE (normal BAM)
```

**Processes:**

| Process | GATK Tool | Purpose |
|---------|-----------|---------|
| `GATKCOV_COUNT` | CollectReadCounts | Read depth in 100 bp bins |
| `GATKCOV_BAF` | CollectAllelicCounts | B-allele frequencies |
| `GATKCOV_CALL` | ModelSegments + CallCopyRatioSegments | Somatic CNV calling |
| `GATKCOV_CALL_GERMLINE` | — | Germline CNV from normal |
| `OVERLAP_GENES` | bedtools | Annotate CNV segments with gene names |
| `FILTER_CNVS_PANEL` | Custom | Keep segments overlapping `params.PANEL_CNV` genes |

**PON files are sex-specific:**
```groovy
params.GATK_PON_FEMALE = "/path/to/female_pon.hdf5"
params.GATK_PON_MALE   = "/path/to/male_pon.hdf5"
```

Sex is inferred from `meta.sex` (set in input CSV or derived from alignment data). `GATKCOV_CALL_GERMLINE` is skipped for tumor-only groups (no normal to call germline segments from) — see [§3.3 Tumor-only Mode](#33-tumor-only-mode).

**Outputs:**
- `tum_plot` — CNV segmentation plot
- `bed` — CNV calls in BED format
- `json` — Structured export for Coyote
- `ascat_segments` — `[group, meta, purity_ploidy.tsv, segments.tsv, segments_raw.tsv]`
- `ascat_loh_heh` — `[group, meta, loh.tsv, subclonal.tsv, gains.tsv, chrom.tsv, heh.txt]`
- `ascat_cn_summary` — `[group, meta, cn_summary.json, arms.tsv, events.tsv]`: the standardised result, meant for reporting. Not yet consumed by `ADD_TO_DB`, `COYOTE_YAML` or `VISUALIZE`.

#### ASCAT: allele-specific copy number, LOH and HeH (paired only)

GATK CNV gives total copy number. ASCAT adds the copies of each parental allele, which is what loss of heterozygosity (LOH) and a karyotype-style chromosome count need. The chain fits purity and ploidy with ASCAT, checks the fit, then summarises it into LOH segments, gains, subclonal events, a chromosome count and a high hyperdiploidy (HeH) verdict, plus a standardised JSON for reporting.

```
tumor+normal CRAM ──► ASCAT_PREPARE ──► ASCAT_FIT ──────────► ASCAT_LOH_HEH
                      alleleCounter     GC/RT correction,      LOH, gains, subclonal events,
                      → LogR / BAF      ASPCF, runAscat,       chromosome count, HeH verdict,
                                        both-allele check, QC  cn_summary.json (schema 1.3)
```

| Process | Script | Purpose | Resources |
|---|---|---|---|
| `ASCAT_PREPARE` | `bin/ascat.R prepare` | alleleCounter at the 1000G SNPs in tumor and normal → LogR/BAF | 24 CPU (one per chromosome), 64 GB; hours |
| `ASCAT_FIT` | `bin/ascat.R fit` (+ `bin/fitlib.R`) | correction, segmentation, purity/ploidy fit, both-allele check, QC of the pair → integer and unrounded allele copies per segment | minutes |
| `ASCAT_LOH_HEH` | `bin/cn_summary.py` (+ the `bin/cnsummary/` package) | summary of the fit (stdlib Python ≥ 3.6) | seconds |

**When it runs**
- Enabled by `params.ascat` (default `true`); container `params.ascat_container`, built from `spike/loh_ascat/ascat.def` (R ≥ 4.3, ASCAT v3 from GitHub, alleleCounter).
- Only groups with a tumor and a normal. Sex comes from the normal's `meta.sex` (`F` → XX, `M` → XY); a group without it is skipped with a warning. Tumor-only groups are skipped (see §3.3).
- Preparation and fit are separate processes, so a refit (other ploidy limits, or `--purity=`/`--ploidy=` in `ASCAT_FIT`'s `ext.args`) reuses the allele counts on `-resume`.

**Code layout.** `bin/ascat.R` holds the ASCAT-specific steps. `bin/fitlib.R` holds what the spike's FACETS script shares with it: argument parsing, the both-allele scoring, `refit_decision`, `prefer_undoubled`, and the arithmetic of the sex check and of tumour cells in the normal. `bin/cn_summary.py` is the command line of the summary. The rules are in `bin/cnsummary/`:

| Module | Contents |
|---|---|
| `genome.py` | GRCh38 chromosomes, centromeres, arms, PARs |
| `tables.py` | TSV and number reading |
| `settings.py` | every threshold: the options and the JSON `params` |
| `segments.py` | segment loading and preparation |
| `calls.py` | LOH, subclonal events, segmental gains |
| `count.py` | chromosome count, range, verdict |
| `layers.py` | arm table, events |
| `qc.py` | purity status, sex, loss of Y, flags |
| `summary.py` | runs the steps in order |
| `report.py` | output files |

The spike (`spike/loh_ascat/`) is the development copy of these files.

##### Prepare (`bin/ascat.R prepare`)

`ascat.prepareHTS` runs alleleCounter at the 1000 Genomes SNPs (`ascat_loci_prefix`, `ascat_alleles_prefix`; ~1 per 1.3 kb) in tumor and normal, one chromosome per thread (1–22, X), with mapping quality ≥ 35 and base quality ≥ 20. It writes `<T>_{tumour,normal}{LogR,BAF}.txt`: LogR is the tumor/normal depth ratio, BAF the alternative-allele fraction. The CRAMs carry no reference path and the cluster has no internet, so the FASTA reaches alleleCounter as `-r` through `additional_allelecounter_flags` (ASCAT 3.2 has no `ref.fasta`); the script stops early if neither argument exists. Sex is passed, not inferred.

##### Fit (`bin/ascat.R fit`)

At SNPs het in the normal, a segment with depth ratio `2^LogR` and major-allele fraction b has, at purity p and tumour ploidy ψ, `M = 2^LogR × (2(1−p) + pψ)` copies in the mixture, and per tumour cell `major = (bM − (1−p))/p`, `minor = ((1−b)M − (1−p))/p`. Fitting means choosing p and ψ.

| # | Step | What it does |
|---|---|---|
| 1 | `ascat.loadData` | reads the four files; het SNPs are those het in the normal |
| 2 | sex check (`sex_check`) | the normal's heterozygosity on X outside the PARs relative to the autosomes: ≥ 0.3 XX, ≤ 0.1 XY, else unclear; compared with the given sex |
| 3 | `ascat.correctLogR` | removes LogR waves that follow GC content (`ascat_gc`) and replication timing (`ascat_rt`); `before_/after_correction_*.png` |
| 4 | `ascat.aspcf` | joint piecewise-constant segmentation of LogR and BAF (BAF at het SNPs, mirrored above 0.5); `<T>.ASPCF.png` |
| 5 | `ascat.runAscat` (γ = 1) | grid over purity 0.10–1.00 and ψ 1.0–5.4; score = distance of each segment's **minor** allele to a whole number, weighted by SNPs, balanced segments at 5%; lowest local minimum with ploidy in `ascat_min_ploidy`–`ascat_max_ploidy` (1.5–4.0); rounds segments to `nMajor`/`nMinor` and keeps the unrounded `nAraw`/`nBraw`; `<T>.sunrise.png`, `<T>.ASCATprofile.png` |
| 6 | refit guard (`aberrant_fraction`) | share of autosomal het SNPs in segments more than 0.15 copies from 1+1 at ASCAT's fit. Below `ascat_min_aberrant` (0.5%, ~15 Mb): `purity_status` = `not identifiable`, no refit, the summary reports events as a share of all cells |
| 7 | both-allele check (`allele_fit`) | same segments and grid, both alleles scored, each allele's squared distance capped at `ascat_allele_cap`² (0.15), negative copies uncapped. A point is feasible when every confident segment (≥ 200 het SNPs, ≥ 1 Mb) has both alleles ≥ −0.1 copies; the purity floor is the lowest feasible purity at ASCAT's ψ. Local minima (7 × 7 window) of the feasible points with ploidy in bounds; `prefer_undoubled` keeps the lower-ploidy minimum unless its doubling is ≥ 1.5× better |
| 8 | refit (`refit_decision`) | `runAscat` rerun at the check's best (`rho_manual`, `psi_manual`) when ASCAT's fit is below the purity floor; or is a doubling (ploidy ≥ 1.6×) of the best without being ≥ 1.5× better; or scores ≥ 1.25× worse, with the best ≥ 1.2× better than the runner-up and purity differing by ≥ 0.05 or ploidy by ≥ 0.2. If the rerun gives no fit, ASCAT's own is kept |
| 9 | tumour in the normal (`tin_estimate`) | the normal's allele shift towards the tumour's major allele in clearly imbalanced segments (a − b ≥ 0.5, ≥ 1 Mb), corrected for SNPs whose tumour allele ratio points the wrong way: the fraction of tumour cells in the normal |
| 10 | `ascat.metrics`, writes | tables below, `<T>.ascat.RData` (`bc`, `out`, `ascat_fit`, `qc`) |

Why steps 6–8 exist:
- **Ploidy cap 4.0**: higher bounds let near-diploid genomes fit a doubled solution. Genomes above ploidy 4 would be forced into a wrong fit.
- **Minor-allele scoring**:
  - A clonal trisomy (2+1) has minor 1 at every purity, and balanced gains count only 5%, so ASCAT's purity is set by the subclonal segments.
  - A large subclonal cnLOH can pull purity down until it reads as clonal, and the trisomies then read as 2.3 copies.
  - Scoring both alleles lets the clonal gains set purity. The cap stops a subclonal segment from costing more the further off it is. Uncapped negative copies let a deletion bound purity from below.
- **Feasibility**: a small clonal deletion needs a high purity, but in the distance it can be outweighed by a large subclonal change. Feasibility enforces it.
- **Doubling**: a doubled solution scores at least as well as its undoubled one, so the lower ploidy is preferred unless the doubling is clearly better.
- **Cap 0.15**: a larger cap lets a balanced gain in most tumour cells pull purity down again.
- **Guard 0.5%**: flat genomes give 0–0.1%; a single clonal 30 Mb deletion gives ~1%.

A manual fit (`--purity=` and `--ploidy=` in `ASCAT_FIT`'s `ext.args`) skips steps 6–8 (`fit_used` = `manual`). When ASCAT finds no fit (non-aberrant samples), purity is `NA`, `fit_used` = `none`, `purity_status` = `not identifiable`, and the process still succeeds.

`purity_ploidy.tsv` holds:
- both fits (`ascat_purity`, `ascat_ploidy`, `allele_purity`, `allele_ploidy`) and `fit_used` (`ascat` / `allele` / `manual` / `none`);
- `purity_floor`, `aberrant_het_fraction`, `purity_status`, and a `note` saying why;
- the QC of the pair: `sex_data`, `sex_het_ratio`, `sex_check`, `tin_fraction`, `tin_se`, `tin_n_het`.

The other tables:
- `segments_raw.tsv`: the unrounded copies per tumor cell (`nAraw`, `nBraw`) and the het SNPs per segment (`n_het`).
- `ownfit.*`: after a refit, ASCAT's own fit is written as well.
- `metrics.tsv`: `ascat.metrics`, i.e. noise (MAPD), GC/RT correction, segment counts, homozygous deletions, LOH fraction, `WGD`, `GI`.

##### Summary (`bin/cn_summary.py`)

ASCAT rounds every segment to whole copies, which hides subclonal events and can give a balanced chromosome an odd state. The summary therefore reads the unrounded segments and works in this order:

| Step | What it does | Why |
|---|---|---|
| Rescale | Divides out the inflation of gained alleles, estimated as the median major allele − 1 of the 2+1 segments (flag above 1.1×) | An under-fitted purity makes clonal trisomies read 2.3 and look subclonal |
| Balanced correction | A segment with equal alleles (\|a − b\| < 0.15) but an odd state becomes even (2+1 → 2+2) | Rounding artefact of a low purity fit |
| Artefact filters | Segments with < 25% of the median het-SNP density (germline homozygosity) and the X-transposed region (chrX:88–93.5 Mb) are left out of LOH and subclonal calls | Missing heterozygosity that the normal shares is not tumour LOH |
| LOH | Segments with minor allele 0, ≥ 1 Mb, typed LOH+loss (1+0), cnLOH (2+0) or LOH+gain (≥ 3+0); neighbours merge only with the same copy number; fraction of tumour cells = 1 − remaining minor allele | A single X in an XY sample is not LOH |
| Subclonal events | Segments whose alleles sit between integers: subclonal cnLOH, LOH or gain, ≥ 5 Mb, with the fraction of tumour cells. Limits: 0.12 below 1 for a lost minor allele, 0.15 for a gain. A gain without allelic shift (k+k, depth only) needs ≥ 20 Mb of such gain on its chromosome | Changes in fewer than ~12–15% of tumour cells are not called; short depth-only gains are noise |
| Segmental gains | Gains with a whole-number state above normal (e.g. 2+1) on ≥ 5 Mb that are neither an arm nor a chromosome gain; subclonal gains inside them are reported there once | A clonal gain of less than half an arm would otherwise not be reported |
| Purity status | `ok`, `not_identifiable`, `lower_bound` (only gains, no loss or LOH: purity is only bounded below), `manual`, `no_fit` | Fractions are of tumour cells when purity is measured, otherwise of all cells |
| Arms | Modal total copy number per arm; arms are split at the hg38 cytoband centromere region, which belongs to neither arm | Fixed rows across samples and tools |
| Chromosome count | Karyotype rule: a chromosome is gained or lost only when both arms share the modal total, each over ≥ 50% of the arm; acrocentrics (13, 14, 15, 21, 22) use the q arm over ≥ 90% (`--whole_frac`). Arm-level changes (del(7q), +1q) are listed, not counted. A chromosome within ±0.05 of its threshold, the count before the balanced correction, or a gain whose counting depends on an unmeasured purity makes the count a range. Y is assumed present for XY | del(7q) must not count as −7 |
| Verdict | `not_applicable` with `--heh_verdict no` (the `solid` profile); `possible_doubled_hypodiploidy` with ≥ 5 whole autosomes in 2+0; `HeH` at 51–67 chromosomes; otherwise `not_HeH` | A doubled near-haploid clone can reach an HeH-range count |
| Events | Chromosome, arm or segment scope; a segment event is dropped when an arm/chromosome event of the same type covers ≥ 90% of it. Clonality: `clonal` ≥ 0.85 of tumour cells, `possibly_subclonal` 0.70–0.85, `subclonal` < 0.70 | Arm fractions scatter by ~0.1, so one cut-off flips labels on noise |

Count classes follow the ALL modal-number groups: near-haploid 24–30, low hypodiploid 31–39, high hypodiploid 40–43, near-diploid hypodiploid 44–45, diploid 46, low hyperdiploid 47–50, high hyperdiploid 51–67, near-triploid 68–80, near-tetraploid > 80.

`ASCAT_LOH_HEH`'s `ext.args` set `--whole_frac 0.9 --masked_hypo_min_loh 5 --min_loh_len 1000000 --subclonal_min_len 5000000 --heh_verdict ${params.ascat_heh_verdict}`. Every other threshold is a default in `bin/cnsummary/settings.py` and is recorded in the JSON `params`.

##### Outputs (`${outdir}/${subdir}/ascat/`)

| File | Content |
|---|---|
| `{id}.ascat.purity_ploidy.tsv` | both fits, `fit_used`, `purity_floor`, `purity_status`, QC of the pair (sex, tumour in normal) |
| `{id}.ascat.segments.tsv`, `.segments_raw.tsv` | integer states; unrounded copies and `n_het` |
| `{id}.ascat.ownfit.*` | after a both-allele refit: ASCAT's own fit (purity_ploidy, segments, segments_raw) |
| `{id}.ascat.heh.txt` | human-readable summary, one item per line (long items continue on indented lines): fit, purity status and floor, sex check, tumour in normal, count (with range), gains and losses, arm-level changes, corrections, segmental gains, LOH, subclonal events, verdict |
| `{id}.ascat.loh.tsv`, `.gains.tsv`, `.subclonal.tsv`, `.chrom.tsv` | LOH segments; segmental gains; subclonal events; per-chromosome mean copy number, LOH fraction and modal state |
| `{id}.ascat.cn_summary.json` | standardised result, schema 1.3. Sections: `sample` (purity only when measured, `purity_status`, `purity_floor`, count and count range, `verdict_code`, sex check, tumour in normal, flags), `arms`, `events` (incl. segmental gains), `artefacts` (germline homozygosity, masked regions, balanced corrections), `params`, `conventions` |
| `{id}.ascat.arms.tsv`, `.events.tsv` | the arm and event layers as tables |
| `{id}.ascat.metrics.tsv` | `ascat.metrics` QC (noise, correction, segments, homozygous deletions, LOH fraction, WGD, GI) |
| `*.png`, `{id}.ascat.RData` | ASCAT plots (profile, ASPCF, sunrise) and R objects |

Field-level definitions of the JSON: `spike/loh_ascat/README.md`. The spike summarises FACETS with the same code (after `to_segments.py`), so both tools share the format.

##### Parameters

| Param | Default | Purpose |
|---|---|---|
| `ascat` | `true` | run the chain |
| `ascat_container` | `${container_dir}/ascat.sif` | ASCAT + alleleCounter image |
| `ascat_loci_prefix`, `ascat_alleles_prefix`, `ascat_gc`, `ascat_rt` | 1000G GRCh38 set | SNP loci/alleles, GC and replication-timing correction (see §8) |
| `ascat_min_ploidy`, `ascat_max_ploidy` | 1.5, 4.0 | ploidy search range |
| `ascat_allele_refit` | `true` | both-allele purity check |
| `ascat_allele_cap` | 0.15 | cap on each allele's distance in that check |
| `ascat_min_aberrant` | 0.005 | refit guard: below this share of het SNPs in aberrant segments, purity is not identifiable |
| `ascat_heh_verdict` | `"yes"`; `"no"` in `solid` | report the HeH verdict (ALL only) |

##### Limitations
- All thresholds are heuristics tuned on the development pairs and checked against PURPLE on a set of AML pairs; none is validated against karyotype or FISH.
- Changes in fewer than ~12–15% of tumour cells are missed.
- A flat genome has no measurable purity (`not_identifiable`). PURPLE gets purity from somatic SNVs there; this chain does not. In the spike, the SNV check `snv_purity.py` can confirm a refit or supply a purity, but it needs a genome-wide somatic VCF, which the pipeline does not publish yet.
- Loss of Y is not assessed in the pipeline (flag `y_not_assessed`). The spike measures it from Y read depth (`y_depth.py`).
- The germline-homozygosity filter judges whole segments, so it cannot trim a low-density stretch inside a long segment.
- Not yet reported to Coyote or Gens: `ascat_cn_summary` is emitted but unused.

**Checks** (in `spike/loh_ascat/`):
- `python3 test_cn_summary.py`: synthetic cases for the summary rules.
- `Rscript test_ascat_qc.R`: base R; the sex check, tumour in the normal, the doubling rule, the refit decision and the both-allele check.
- `regress.py` and `regress_fit.R`: compare every result on a collect before and after a code change.

---

### 6.8 SV_CALLING

**File:** `subworkflows/local/sv_calling.nf`

Structural variant detection with Manta, followed by fusion gene annotation.

**Process chain:**

```
bam (tumour+normal) ──► MANTA ──► MANTA_SV ──► SNPEFF ──► FILTER_FUSIONS_PANEL ──► COMBINE_FUSIONS
                                              └──► SNPEFF_SV_ANN
```

**Processes:**

| Process | Tool | Purpose |
|---------|------|---------|
| `MANTA` | Manta | Genome-wide SV calling (BND, DEL, DUP, INV, INS) |
| `MANTA_SV` | Custom | Extract specific SV subtypes |
| `SNPEFF` | SnpEff | Annotate predicted fusion genes |
| `SNPEFF_SV_ANN` | SnpEff | Annotate remaining SV types |
| `FILTER_FUSIONS_PANEL` | Custom | Keep fusions in `params.PANEL_FUS` gene list |
| `COMBINE_FUSIONS` | Custom | Merge all SV types into final output |

**Outputs:**
- `fusions` — Final fusion and SV calls (VCF + TSV)

For tumor-only groups, `MANTA` runs in `--tumorBam`-only mode and emits `tumorSV.vcf.gz`; the `manta_vcf_normal` output is `optional: true` and absent — see [§3.3 Tumor-only Mode](#33-tumor-only-mode).

---

### 6.9 VISUALIZE

**File:** `subworkflows/local/visualize.nf`

Generates visualisations from CNV and germline variant data. Integrates with `GENS` for genomic coverage plotting. Triggered after CNV calling and DNAscope germline calling are complete.

---

### 6.10 ADD_TO_DB

**File:** `subworkflows/local/add_to_db.nf`

Loads final variant calls into the Coyote clinical database.

**Processes:**

| Process | Purpose |
|---------|---------|
| `COYOTE` | Import SNV, CNV, SV VCFs to Coyote |
| `COYOTE_YAML` | Export structured YAML for database ingestion |

**Control parameters:**

| Param | Default | Description |
|-------|---------|-------------|
| `params.coyote_group` | `"tumwgs"` | Coyote database group |
| `params.assay` | `"tumwgs"` | Assay identifier |
| `params.cdm` | `"tumwgs"` | CDM identifier |

Set `params.coyote = false` to skip database import during testing.

---

## 7. Profiles

### `hema` — Hematologic malignancies

```groovy
params.PANEL_SNV  = "20250321_Hema_snv_genes"
params.PANEL_CNV  = "20250917_Hema_cnv_addon_genes"
params.PANEL_FUS  = "20250919_Hema_fusion_genes"
params.coyote_group = "tumwgs-hema"
```

### `solid` — Solid tumours

```groovy
params.PANEL_SNV  = "20250610_BTB_snv.panel"
params.PANEL_CNV  = "20250610_BTB_cna.panel"
params.PANEL_FUS  = "20250610_BTB_fusion.panel"
params.coyote_group = "tumwgs-solid"
```

### `test` — Development

Minimal gene panels, reduced resource allocations. Suitable for testing on a subset of data or in a local environment.

### `hopper` / `trannel` — Cluster environments

Sets SLURM executor, queue, Singularity container path, and resource limits for the Grace/Hopper or Trannel compute clusters.

---

## 8. Reference Data & Annotation Databases

All reference data paths are set in `nextflow.config` under `params`. When updating any database version, update both the path and the version table in this document.

| Resource | Parameter | Version |
|----------|-----------|---------|
| Reference genome | `params.genome_file` | GRCh38 GCA_000001405.15 (no alt) |
| Known indels | `params.known` | Mills + 1000G gold standard |
| dbSNP | `params.dbSnp` | v146 |
| VEP cache | `params.VEP_CACHE` | v113.0 |
| gnomAD exomes | `params.GNOMAD` | v4.0 |
| CADD scores | `params.CADD` | v1.7 |
| COSMIC | `params.COSMIC` | v92 |
| GATK genomic bins | `params.GATK_intervals_full` | 100 bp bins |
| PON SNV (Freebayes/VarDict) | `params.PON` | 206 normals |
| PON CNV (female) | `params.GATK_PON_FEMALE` | — |
| PON CNV (male) | `params.GATK_PON_MALE` | — |
| Gene annotations | `params.GENE_BED` | Gencode v33 protein-coding |
| Exon coordinates | `params.EXON_BED` | Ensembl 98 + ClinVar 5 bp pad |
| ASCAT loci / alleles | `params.ascat_loci_prefix`, `params.ascat_alleles_prefix` | 1000G GRCh38 (sarek/igenomes set), `chr` stripped |
| ASCAT GC / replication timing | `params.ascat_gc`, `params.ascat_rt` | 1000G GRCh38 (sarek/igenomes set), `chr` stripped |

---

## 9. Software Versions

| Tool | Version | Stage |
|------|---------|-------|
| Sentieon | 202308.03 / 202503 | Alignment, dedup, BQSR, QC, TNscope, DNAscope |
| GATK | 4.1.9.0 / 4.2.x | CNV calling |
| ASCAT + alleleCounter | v3 (GitHub) | Allele-specific CN, LOH, HeH (`ascat.sif`) |
| Ensembl VEP | 113.0 | SNV annotation |
| DeepSomatic | 1.9.0 | SNV calling |
| Manta | latest stable | SV calling |
| Pindel | latest stable | Indel calling |
| FreeBayes | latest stable | SNV calling |
| VarDict | latest stable | SNV calling |
| SnpEff | latest stable | SV annotation |
| fastp | latest stable | FASTQ trimming |
| seqtk | latest stable | Downsampling |
| bcftools / vcftools | latest stable | VCF processing |
| SVDB | latest stable | SV database merging |
| bedtools | latest stable | Region operations |

All tools run inside a Singularity container (`tumwgs_container.sif`). Software versions are automatically collected at the end of each run and written to `software_versions.yml`.

---

## 10. Output Structure

```
{params.outdir}/tumwgs/
├── bam/
│   ├── {id}.bam                    # Analysis-ready BAM (BQSR)
│   └── {id}.bam.bai
├── cram/
│   ├── {id}.dedup.cram             # Deduplicated CRAM (archive)
│   └── {id}.bqsr.cram              # BQSR CRAM
├── vcf/
│   ├── {id}.agg.vcf                # Aggregated multi-caller VCF
│   ├── {id}.vep.vcf                # VEP-annotated VCF
│   ├── {id}.filtered.vcf           # Final filtered VCF (panel)
│   ├── {id}.dnascope.vcf           # Germline variants
│   └── {id}.tnscope.vcf.gz         # Raw TNscope calls
├── cnv/
│   ├── {id}.cnv.bed                # CNV segments (BED)
│   ├── {id}.cnv.vcf                # CNV in VCF format
│   └── {id}.cnv.json               # Structured CNV for Coyote
├── ascat/                          # paired groups only; {id} = tumor id
│   ├── {id}.ascat.purity_ploidy.tsv
│   ├── {id}.ascat.segments.tsv     # nMajor/nMinor per segment (segments_raw.tsv: unrounded)
│   ├── {id}.ascat.ownfit.*         # after a both-allele refit: ASCAT's own fit
│   ├── {id}.ascat.loh.tsv          # LOH segments >= 1 Mb, typed, with tumor-cell fraction
│   ├── {id}.ascat.gains.tsv        # segmental gains >= 5 Mb (part of an arm)
│   ├── {id}.ascat.subclonal.tsv    # subclonal cnLOH / LOH / gains >= 5 Mb
│   ├── {id}.ascat.chrom.tsv        # per-chromosome mean CN, LOH fraction, modal state
│   ├── {id}.ascat.heh.txt          # chromosome count + HeH verdict
│   ├── {id}.ascat.cn_summary.json  # standardised summary (schema 1.3): sample, arms, events, artefacts, params
│   ├── {id}.ascat.arms.tsv         # fixed per-arm table (1p..Xq): state, class, LOH fraction
│   ├── {id}.ascat.events.tsv       # reportable events: chromosome / arm / segment scope, clonality, cell fractions
│   └── *.png, {id}.ascat.RData     # ASCAT plots (incl. sunrise), R objects
├── sv/
│   ├── {id}.fusions.vcf            # Fusion gene calls
│   └── {id}.sv.annotated.vcf       # All annotated SVs
├── QC/
│   ├── {id}.qc.txt                 # Alignment QC summary
│   ├── {id}.dedup.metrics          # Duplication rate
│   └── {id}.idsnp.txt              # Identity SNP results
├── plots/
│   └── {id}.cnv_plot.png           # CNV segmentation plot
├── {id}.coyote.yml                 # Coyote import YAML
└── software_versions.yml           # All tool versions
```

---

## 11. Adding New Components

### Adding a New SNV Caller

1. Create `modules/local/<caller>/main.nf` following the module template (see `CLAUDE.md`).
2. Add resource config in `configs/modules/snv_calling.config`.
3. Import and invoke in `subworkflows/local/snv_calling.nf`.
4. Feed output VCF into `CONCATENATE_VCFS` → `AGGREGATE_VCFS`.
5. Add a toggle param in `nextflow.config`: `params.<caller> = true`.
6. Add version collection to `CUSTOM_DUMPSOFTWAREVERSIONS`.
7. Update §9 (Software Versions) in this document.

### Adding a New Gene Panel

1. Place the panel file in the appropriate reference data location.
2. Add `params.PANEL_XXX = "/path/to/panel"` in `nextflow.config` under the target profile.
3. Pass the param to the relevant filter process in the subworkflow.

### Adding a New Profile

1. Add a `profiles { newprofile { ... } }` block in `nextflow.config`.
2. At minimum set: `params.PANEL_SNV`, `params.PANEL_CNV`, `params.PANEL_FUS`, `params.coyote_group`.
3. Document the new profile in §7.

### Updating a Reference Database

1. Update the file path in `nextflow.config` (`params.*`).
2. Update the version in §8 of this document.
3. Note the change in `CHANGELOG.md`.

---

## 12. Parameter Reference

### Core Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `params.csv` | — | Input CSV file path (required) |
| `params.outdir` | — | Output directory |
| `params.subdir` | `""` | Subdirectory within outdir |
| `params.genome_file` | — | Reference FASTA path |
| `params.crondir` | — | Cron log directory |

### Alignment

| Parameter | Default | Description |
|-----------|---------|-------------|
| `params.bwa_shards` | `8` | Parallel BWA shards |
| `params.K_size` | `100000000` | BWA reads per block |
| `params.known` | — | Known indels VCF |
| `params.dbSnp` | — | dbSNP VCF |

### SNV Calling

| Parameter | Default | Description |
|-----------|---------|-------------|
| `params.freebayes` | `true` | Enable Freebayes |
| `params.vardict` | `true` | Enable VarDict |
| `params.tnscope` | `true` | Enable Sentieon TNscope |
| `params.deepsomatic` | `true` | Enable DeepSomatic |
| `params.pindel` | `true` | Enable Pindel |
| `params.dnascope` | `true` | Enable germline calling |
| `params.fb_var_freq_cutoff_p` | `'0.03'` | Freebayes min VAF (paired) |
| `params.vardict_var_freq_cutoff_p` | `'0.03'` | VarDict min VAF (paired) |
| `params.tnscope_var_freq_cutoff_p` | `'0.01'` | TNscope min VAF (paired) |

### SNV Annotation & Filtering

| Parameter | Default | Description |
|-----------|---------|-------------|
| `params.PANEL_SNV` | profile-specific | SNV gene panel file |
| `params.VEP_CACHE` | — | VEP cache directory |
| `params.CADD` | — | CADD score file |
| `params.GNOMAD` | — | gnomAD frequency file |
| `params.COSMIC` | — | COSMIC mutations file |
| `params.filter_freq` | `0.05` | Max gnomAD AF |
| `params.filter_field_filter` | `"FAIL_PON*,..."` | FILTER tags to exclude |
| `params.override_filter_terms` | `"CLIN_SIG=likely_pathogenic,pathogenic"` | Filter override |

### CNV Calling

| Parameter | Default | Description |
|-----------|---------|-------------|
| `params.gatk_cnv` | `true` | Enable GATK CNV |
| `params.GATK_PON_FEMALE` | — | Female CNV PON (HDF5) |
| `params.GATK_PON_MALE` | — | Male CNV PON (HDF5) |
| `params.GATK_intervals_full` | — | 100 bp genomic bins |
| `params.PANEL_CNV` | profile-specific | CNV gene panel file |

### SV Calling

| Parameter | Default | Description |
|-----------|---------|-------------|
| `params.PANEL_FUS` | profile-specific | Fusion gene panel file |
| `params.FUSIONS_CNV` | — | Fusion-associated SV definitions |

### Database & Output

| Parameter | Default | Description |
|-----------|---------|-------------|
| `params.coyote_group` | `"tumwgs"` | Coyote DB group |
| `params.assay` | `"tumwgs"` | Assay name |
| `params.cdm` | `"tumwgs"` | CDM identifier |
| `params.coyote` | `true` | Enable DB import |

### Sampling

| Parameter | Default | Description |
|-----------|---------|-------------|
| `params.sample` | `true` | Enable downsampling |
| `params.trimfq` | `false` | Enable fastp trimming |
