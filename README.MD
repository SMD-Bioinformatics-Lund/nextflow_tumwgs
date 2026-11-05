# SMD-Bioinformatics Lund – Somatic Whole Genome Sequencing (SWGS) Pipeline

This page describes how to run the Somatic Whole Genome Sequencing Pipeline (main.nf) on the SMD-Bioinformatics computing cluster (Hopper/Grace), as well as how the pipeline integrates with the automated infrastructure.

## Overview

The main entry point for the pipeline is:

```
main.nf
```

Before running the pipeline, make sure to update the Nextflow config (nextflow.config) to match the server environment and to choose the correct panel based on the biological question—e.g., hematology or solid tumor.

Our compute environments currently include the Hopper and Grace SLURM clusters.
A corresponding Nextflow profile (hopper) is available in the configuration file.

The elaborate technical documentation on the bioinformatics softwares used are described in the [documenation](./doc).

### Pipeline Overview

```mermaid
flowchart TD
    csv[/"Samplesheet CSV<br/>fastq / bam / cram"/] --> CHECK_INPUT

    subgraph PRE["Pre-processing"]
        direction TB
        CHECK_INPUT --> SAMPLE["SAMPLE<br/>seqtk · fastp"]
        SAMPLE --> ALIGN_SENTIEON["ALIGN_SENTIEON<br/>BWA (8 shards) · Dedup · BQSR"]
        CHECK_INPUT --> ALIGN_BAM_CRAM["ALIGN_BAM_CRAM<br/>bam_state: raw / processed"]
    end

    ALIGN_SENTIEON --> BAMS(["BQSR bam · dedup cram"])
    ALIGN_BAM_CRAM --> BAMS

    BAMS --> BAM_QC["BAM_QC<br/>Sentieon QC · Somalier · VerifyBamID2 · ID-SNP"]
    BAMS --> DUX4["DUX4IGH_CALLING<br/>Pelops"]
    BAMS --> SNV_CALLING["SNV_CALLING<br/>TNscope ML · DeepSomatic · FreeBayes<br/>VarDict · Pindel · DNAscope"]
    BAMS --> CNV_CALLING["CNV_CALLING<br/>GATK CNV · BAF"]
    BAMS --> SV_CALLING["SV_CALLING<br/>Manta · snpEff · fusion panels"]

    SNV_CALLING -- "DNAscope GVCF → GVCFtyper" --> GERMLINE["GERMLINE_ANNOTATE<br/>VEP + ClinVar · mark_germlines"]
    SNV_CALLING -- "aggregated VCF" --> SNV_ANNOTATE["SNV_ANNOTATE<br/>PON · VEP 113 · panel filter"]
    GERMLINE --> SNV_ANNOTATE

    SNV_CALLING --> VISUALIZE["VISUALIZE<br/>GENS v4"]
    CNV_CALLING --> VISUALIZE

    BAM_QC --> CDM[("CDM<br/>QC metrics")]
    SNV_ANNOTATE --> ADD_TO_DB
    CNV_CALLING --> ADD_TO_DB
    SV_CALLING --> ADD_TO_DB
    ADD_TO_DB["ADD_TO_DB"] --> COYOTE[("Coyote3")]
    VISUALIZE --> GENS[("GENS")]

    classDef input fill:#24B064,stroke:#1a7a46,color:#fff
    classDef proc fill:#f2f9f5,stroke:#24B064,color:#000
    classDef db fill:#0DC09D,stroke:#087a63,color:#fff
    class csv,BAMS input
    class CHECK_INPUT,SAMPLE,ALIGN_SENTIEON,ALIGN_BAM_CRAM,BAM_QC,DUX4,SNV_CALLING,CNV_CALLING,SV_CALLING,GERMLINE,SNV_ANNOTATE,VISUALIZE,ADD_TO_DB proc
    class CDM,COYOTE,GENS db
```

### Running the Pipeline Manually

#### 1. Hematology Pipeline

```bash
module load singularity Java nextflow/21.10.6

nextflow run main.nf \
    -entry SWGP \
    -c nextflow.config \
    --csv sample.csv \
    -profile hema \
    -with-report /path/to/workdir/nextflow/reports/sample.report.html \
    -with-trace /path/to/workdir/nextflow/reports/sample.trace.txt \
    -with-timeline /path/to/workdir/nextflow/reports/sample.timeline.html \
    -work-dir /path/to/workdir/nextflow_tmp
```

#### 2. Solid Tumor Pipeline

```bash
module load singularity Java nextflow/21.10.6

nextflow run main.nf \
    -entry SWGP \
    -c nextflow.config \
    --csv sample.csv \
    -profile solid \
    -with-report /path/to/workdir/nextflow/reports/sample.report.html \
    -with-trace /path/to/workdir/nextflow/reports/sample.trace.txt \
    -with-timeline /path/to/workdir/nextflow/reports/sample.timeline.html \
    -work-dir /path/to/workdir/nextflow_tmp
```

#### 3. Test / Local Minimal Profile

A lightweight test profile exists for local or simplified runs:

```bash
module load singularity Java nextflow/21.10.6

nextflow run main.nf \
    -entry SWGP \
    -c nextflow.config \
    -profile hopper \
    --csv sample.csv \
    -profile test \
    --dev \
    -with-report /path/to/workdir/nextflow/reports/sample.report.html \
    -with-trace /path/to/workdir/nextflow/reports/sample.trace.txt \
    -with-timeline /path/to/workdir/nextflow/reports/sample.timeline.html \
    -work-dir /path/to/workdir/nextflow_tmp
```

### Integration With SMD Automated Pipeline System

In production, pipeline execution is automated using the
bnf-infrastructure
repository.

Two key components handle this:

```
1. start_nextflow_analysis.pl
```

A Perl script that monitors and triggers pipeline runs based on the availability of CSV files from earlier `Bjorn` system

```
2. pipeline_files.config
```

Defines which pipeline, container, and server settings to use.

Example configuration:

```
[tumwgs-hema]
pipeline = /production_pipeline/nextflow_tumwgs/main.nf -entry SWGP --profile hema
container = /production_pipeline/nextflow_tumwgs/container/tumwgs_container.sif
singularity_version = 3.8.0
nextflow_version = 21.04.2
executor = slurm
cluster = grace
queue = normal

[tumwgs-solid]
pipeline = /production_pipeline/nextflow_tumwgs/main.nf -entry SWGP --profile solid
container = /production_pipeline/nextflow_tumwgs/container/tumwgs_container.sif
singularity_version = 3.8.0
nextflow_version = 21.04.2
executor = slurm
cluster = grace
queue = normal

```
