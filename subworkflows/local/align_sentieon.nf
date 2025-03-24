#!/usr/bin/env nextflow

include { BWA_ALIGN_SHARD               } from '../../modules/local/sentieon/main'
include { BWA_MERGE_SHARDS              } from '../../modules/local/sentieon/main'
include { BAM_CRAM                      } from '../../modules/local/sentieon/main'
include { MARKDUP                       } from '../../modules/local/sentieon/main'
include { REALIGN_INDEL_BQSR            } from '../../modules/local/sentieon/main'
include { CRAM_TO_BAM                   } from '../../modules/local/sentieon/main'

workflow ALIGN_SENTIEON {
    take: 
        bwa_shards          // channel : [0,1...7]          
        fastq_input         // channel: [mandatory] [ val(shard), val(meta), [ reads ] ]
        meta                // channel: [mandatory] [ [sample_id, group, sex, phenotype, paternal_id, maternal_id, case_id] ]

    main:
        ch_versions = Channel.empty()

        BWA_ALIGN_SHARD ( bwa_shards.combine(fastq_input) )
        ch_versions = ch_versions.mix(BWA_ALIGN_SHARD.out.versions)
        
        BWA_MERGE_SHARDS ( BWA_ALIGN_SHARD.out.shard_bam.groupTuple( by: [0,1,2]) )
        ch_versions = ch_versions.mix(BWA_MERGE_SHARDS.out.versions)
        
        BAM_CRAM (  BWA_MERGE_SHARDS.out.merged_bam )
        ch_versions = ch_versions.mix(BAM_CRAM.out.versions)

        MARKDUP ( BAM_CRAM.out.cram_merged  )
        ch_versions = ch_versions.mix(BAM_CRAM.out.versions)

        REALIGN_INDEL_BQSR (    MARKDUP.out.cram_bqsr    )
        ch_versions = ch_versions.mix(REALIGN_INDEL_BQSR.out.versions)

        CRAM_TO_BAM ( REALIGN_INDEL_BQSR.out.cram_bqsr )
        ch_versions = ch_versions.mix(CRAM_TO_BAM.out.versions)
        
    emit:
        bam_bqsr        = CRAM_TO_BAM.out.bam_bqsr                    // channel: [ val(group), val(meta), file(bam), file(bai) ]
        cram_bqsr       = REALIGN_INDEL_BQSR.out.cram_varcall         // channel: [ val(group), val(meta), file(bam), file(bai), file(bqsr.table) ]
        cram_dedup      = REALIGN_INDEL_BQSR.out.cram_bqsr            // channel: [ val(group), val(meta), file(bam), file(bai)] 
        dedup_metrics   = MARKDUP.out.cram_metric
        versions        = ch_versions             // channel: [ file(versions) ]
}