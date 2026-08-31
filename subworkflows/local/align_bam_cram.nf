#!/usr/bin/env nextflow

include { BAM_CRAM                       } from '../../modules/local/sentieon/main'
include { CRAM_INDEX                     } from '../../modules/local/sentieon/main'
include { MARKDUP                        } from '../../modules/local/sentieon/main'
include { REALIGN_INDEL_BQSR             } from '../../modules/local/sentieon/main'
include { CRAM_TO_BAM                    } from '../../modules/local/sentieon/main'
include { CRAM_TO_BAM as CRAM_TO_BAM_PROCESSED } from '../../modules/local/sentieon/main'
include { QUALCAL_ONLY                   } from '../../modules/local/sentieon/main'
include { DEDUP_METRICS_PLACEHOLDER      } from '../../modules/local/sentieon/main'

// Lets a samplesheet row point directly at a bam or cram instead of fastq, so
// already-aligned data can skip BWA/BWA_MERGE_SHARDS. meta.bam_state (set in
// create_meta.nf, default 'raw') decides how much of the usual post-alignment
// chain still needs to run for that row:
//   'raw'       - alignment-only bam/cram: still needs Dedup + Realign + BQSR,
//                 i.e. the same tail ALIGN_SENTIEON runs after BWA.
//   'processed' - already dedup+realign+BQSR'd (e.g. reused output from a
//                 previous run of this pipeline): only the small pieces the
//                 rest of the pipeline actually needs are synthesized (a BQSR
//                 table for the variant callers, a placeholder dedup-metrics
//                 file for BAM_QC's join).
workflow ALIGN_BAM_CRAM {
    take:
        bam     // channel: [ val(group), val(meta), file(bam), file(bai) ]
        cram    // channel: [ val(group), val(meta), file(cram), file(crai) ]

    main:
        ch_versions = Channel.empty()

        // Normalize both input shapes to a uniform [group, meta, cram, crai, bai]
        // before branching on bam_state, so every process below runs exactly once.
        BAM_CRAM ( bam )
        ch_versions = ch_versions.mix(BAM_CRAM.out.versions)

        CRAM_INDEX ( cram )
        ch_versions = ch_versions.mix(CRAM_INDEX.out.versions)

        cram_normalized = BAM_CRAM.out.cram_merged.mix( CRAM_INDEX.out.cram_bai )

        cram_raw       = cram_normalized.filter { g, m, c, i, b -> m.bam_state != 'processed' }
        cram_processed = cram_normalized.filter { g, m, c, i, b -> m.bam_state == 'processed' }

        // --- raw: Dedup + Realign + BQSR, mirroring ALIGN_SENTIEON's tail ---
        MARKDUP ( cram_raw )
        ch_versions = ch_versions.mix(MARKDUP.out.versions)

        REALIGN_INDEL_BQSR ( MARKDUP.out.cram_bqsr )
        ch_versions = ch_versions.mix(REALIGN_INDEL_BQSR.out.versions)

        CRAM_TO_BAM ( REALIGN_INDEL_BQSR.out.cram_bqsr )
        ch_versions = ch_versions.mix(CRAM_TO_BAM.out.versions)

        // --- processed: already dedup+realign+BQSR'd, just fill in the gaps ---
        QUALCAL_ONLY ( cram_processed )
        ch_versions = ch_versions.mix(QUALCAL_ONLY.out.versions)

        DEDUP_METRICS_PLACEHOLDER ( cram_processed.map { g, m, c, i, b -> tuple(g, m) } )
        ch_versions = ch_versions.mix(DEDUP_METRICS_PLACEHOLDER.out.versions)

        CRAM_TO_BAM_PROCESSED ( cram_processed )
        ch_versions = ch_versions.mix(CRAM_TO_BAM_PROCESSED.out.versions)

    emit:
        bam_bqsr        = CRAM_TO_BAM.out.bam_bqsr.mix( CRAM_TO_BAM_PROCESSED.out.bam_bqsr )
        cram_dedup      = REALIGN_INDEL_BQSR.out.cram_bqsr.mix( cram_processed )
        cram_bqsr       = REALIGN_INDEL_BQSR.out.cram_varcall.mix( QUALCAL_ONLY.out.cram_varcall )
        dedup_metrics   = MARKDUP.out.cram_metric.mix( DEDUP_METRICS_PLACEHOLDER.out.cram_metric )
        versions        = ch_versions
}
