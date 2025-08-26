#!/usr/bin/env nextflow

include { SENTIEON_QC          } from '../../modules/local/sentieon/main'
include { COLLECT_QC           } from '../../modules/local/sentieon/main'
include { QC_TO_CDM            } from '../../modules/local/qc/main'
include { ALLELE_CALL          } from '../../modules/local/idSnp/main'
include { SNP_CHECK            } from '../../modules/local/idSnp/main'
include { PAIRGEN_CDM          } from '../../modules/local/idSnp/main'

workflow BAM_QC {
    take:
        bam_bqsr        // channel: [ val(group), val(meta), file(bam), file(bai) ]
        bam_dedup       // channel: [ val(group), val(meta), file(cram), file(crai), file(bai) ]
        dedup_metrics   // channel: [ val(group), val(meta), file(dedup_metrics) ]

    main:
        ch_versions = Channel.empty()
        
        SENTIEON_QC ( bam_dedup.join(dedup_metrics, by:[0,1]) )
        ch_versions = ch_versions.mix(SENTIEON_QC.out.versions)

        COLLECT_QC ( SENTIEON_QC.out.qc )
        ch_versions = ch_versions.mix(COLLECT_QC.out.versions)

        QC_TO_CDM ( COLLECT_QC.out.qc_cdm )


        // Check genotypes of ID-SNPs
        ALLELE_CALL (bam_bqsr)
        ch_versions = ch_versions.mix(ALLELE_CALL.out.versions)

        SNP_CHECK(ALLELE_CALL.out.sample_id_genotypes.groupTuple())
        ch_versions = ch_versions.mix(SNP_CHECK.out.versions)

        PAIRGEN_CDM (SNP_CHECK.out.idsnp_checked)

    emit:
        qcdone                  =   QC_TO_CDM.out.cdm_done    // channel: [ val(group), val(meta), file
        versions                =   ch_versions              // channel: [ file(versions) ]
        dedup_bam_is_metrics    =   COLLECT_QC.out.qc_cdm   // channel: [ val(group), val(meta), file(is_metrics.txt) ] 
}