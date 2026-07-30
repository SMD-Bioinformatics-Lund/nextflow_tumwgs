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

---

## 9. Software Versions

| Tool | Version | Stage |
|------|---------|-------|
| Sentieon | 202308.03 / 202503 | Alignment, dedup, BQSR, QC, TNscope, DNAscope |
| GATK | 4.1.9.0 / 4.2.x | CNV calling |
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
