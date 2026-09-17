#!/usr/bin/env nextflow

include { PON_FILTER               } from '../../modules/local/filters/main'
include { ANNOTATE_VEP             } from '../../modules/local/filters/main'
include { FILTER_PANEL             } from '../../modules/local/filters/main'
include { FIX_VEP                   } from '../../modules/local/filters/main'
include { POST_ANNOTATION_FILTERS  } from '../../modules/local/filters/main'

workflow SNV_ANNOTATE {
    take:
        agg_vcf         // channel: [mandatory] [ val(group), val(meta), file(agg.vcf) ]
        concat_vcfs     // channel: [mandatory] [ val(group), val(vc), file(vcf.gz) ]
        meta            // channel: [mandatory] [ [sample_id, group, sex, phenotype, paternal_id, maternal_id, case_id] ]

    main:
        ch_versions = Channel.empty()

        // Filter with PoN, annotate with VEP, mark germlines
        PON_FILTER { agg_vcf }
        ch_versions = ch_versions.mix(PON_FILTER.out.versions)

        ANNOTATE_VEP { PON_FILTER.out.vcf_pon }
        ch_versions = ch_versions.mix(ANNOTATE_VEP.out.versions)

        FILTER_PANEL { ANNOTATE_VEP.out.vcf_vep }
        ch_versions = ch_versions.mix(FILTER_PANEL.out.versions)

        FIX_VEP { FILTER_PANEL.out.vcf_panel }
        ch_versions = ch_versions.mix(FIX_VEP.out.versions)

        // add filters and mark duplicates
        POST_ANNOTATION_FILTERS { FIX_VEP.out.fixed_vcf }
        ch_versions = ch_versions.mix(POST_ANNOTATION_FILTERS.out.versions)

    emit:
        annotated_variants  =   FILTER_PANEL.out.vcf_panel                      // channel: [ val(group), val(vc), file(vcf.gz) ]
        finished_vcf        =   POST_ANNOTATION_FILTERS.out.filtered_vcf        // channel: [ val(group), val(vc), file(vcf.gz) ]
        versions            =   ch_versions                                     // channel: [ file(versions) ]
}
