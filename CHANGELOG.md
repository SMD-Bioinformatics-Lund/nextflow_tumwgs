# CHANGELOG

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

### v3.0.0

Changes since v2.0_production.

#### Features

    *   DSL2 restructure aligned with Somatic Panel Pipeline; single `SWGP` entry point
    *   Separate `hema` / `solid` profiles; configs split per subworkflow
    *   New SNV callers: TNscope (+ML model), DeepSomatic, Pindel
    *   Germline reporting: DNAscope -> GVCFtyper -> VEP (+ClinVar) -> mark_germlines, merged into somatic VCF for CEBPA 
    *   BAM/CRAM input (skip alignment) via `bam_state` column
    *   DUX4-IGH rearrangement caller (Pelops)
    *   SV/tandem-dup based fusion detection; updated SMD/sarcoma fusion panels
    *   QC: Somalier (sex, pairing), contamination metrics, QC export to CDM
    *   Coyote3 import (CNV json, clarity id) and GENS v4
    *   Software versions collected for all processes
    *   Clinically validated with the tumwgs samples that was earlier run in the research and development project

#### Updates

    *   VEP 113; refreshed annotation databases and curated SNV/CNV gene lists
    *   Masked hg38 genome build
    *   GATK CNV tuned for NovaSeq X
    *   Freebayes/VarDict called on dedup BAM; filters tuned for low-VAF SNVs

#### Fixes

    *   TNscope sharding removed (broke downstream steps)
    *   Tumor-only fixes for Manta, Somalier, TNscope filters
    *   Pindel tumor/normal index swap; ID-SNP pairing; QC subassay field
    *   Coyote output/GENS link fixes

### Release date

    * version 2.0.0 on September 30 2022

#### Features

    *	DSL2 based implemetation of tumor WGS
    *	Able to adabpt to both Hematology and Solid WGS questions
    * 	GENS for both normal and tumor
    * 	Able to handle both tumor/normal and tumor only analysis
    * 	COYOTE can handle unfiltered SNV list but tagged SNV from gene panel can be selected

#### Fixes

    * 	Sharding from 32 to 8 for better I/O specially with bqsr steps in senteion
    * 	Changes with the profile configs
    * 	Changes with the coyote to include sample-wgs
