#!/usr/bin/env nextflow

include { GATKCOV_BAF                          } from '../../modules/local/GATK/main'
include { GATKCOV_COUNT                        } from '../../modules/local/GATK/main'
include { GATKCOV_CALL                         } from '../../modules/local/GATK/main'
include { GATKCOV_CALL_GERMLINE                } from '../../modules/local/GATK/main'
include { OVERLAP_GENES                        } from '../../modules/local/filters/main'
include { FILTER_CNVS_PANEL                    } from '../../modules/local/filters/main'


workflow CNV_CALLING {
    take: 
        cram_dedup            // channel: [mandatory] [ val(group), val(meta), file(cram), file(crai), file(bai) ]
        meta                 // channel: [mandatory] [ [sample_id, group, sex, phenotype, paternal_id, maternal_id, case_id] ]

    main:
        ch_versions = Channel.empty()

        //////////////////////////// GATK SEGMENT CALLING /////////////////////////////////////
        // Do calling somatic CNV calling, use normal allelic counts for somatic as well     //
        ///////////////////////////////////////////////////////////////////////////////////////
        GATKCOV_BAF (  cram_dedup  )
        ch_versions = ch_versions.mix(GATKCOV_BAF.out.versions)

        GATKCOV_COUNT ( cram_dedup )
        ch_versions = ch_versions.mix(GATKCOV_COUNT.out.versions)

        GATKCOV_CALL { GATKCOV_BAF.out.gatk_baf.join(GATKCOV_COUNT.out.gatk_count,by:[0,1]).groupTuple() }
        ch_versions = ch_versions.mix(GATKCOV_CALL.out.versions)

        GATKCOV_CALL_GERMLINE { GATKCOV_BAF.out.gatk_baf.join(GATKCOV_COUNT.out.gatk_count,by:[0,1]).groupTuple() }
        ch_versions = ch_versions.mix(GATKCOV_CALL.out.versions)

        OVERLAP_GENES { GATKCOV_CALL.out.gatkcov_called }
        ch_versions = ch_versions.mix(OVERLAP_GENES.out.versions)

        FILTER_CNVS_PANEL { OVERLAP_GENES.out.annotated_bed }
        ch_versions = ch_versions.mix(OVERLAP_GENES.out.versions)

    emit:
        tum_plot    =   GATKCOV_CALL.out.gatkcov_plot
        norm_plot   =   GATKCOV_CALL_GERMLINE.out.gatkcov_germlline_plot
        count       =   GATKCOV_COUNT.out.gatk_count
        bed         =   FILTER_CNVS_PANEL.out.vcf_panel_bed
        json        =   FILTER_CNVS_PANEL.out.vcf_panel_json
        versions    =   ch_versions                     // channel: [ file(versions) ]

}