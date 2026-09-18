#!/usr/bin/env nextflow

include { PON_FILTER               } from '../../modules/local/filters/main'
include { ANNOTATE_VEP             } from '../../modules/local/filters/main'
include { FILTER_PANEL             } from '../../modules/local/filters/main'
include { FIX_VEP                   } from '../../modules/local/filters/main'
include { POST_ANNOTATION_FILTERS  } from '../../modules/local/filters/main'
include { COMBINE_SOMATIC_GERMLINE } from '../../modules/local/filters/main'

workflow SNV_ANNOTATE {
    take:
        agg_vcf             // channel: [mandatory] [ val(group), val(meta), file(agg.vcf) ]
        concat_vcfs         // channel: [mandatory] [ val(group), val(vc), file(vcf.gz) ]
        meta                // channel: [mandatory] [ [sample_id, group, sex, phenotype, paternal_id, maternal_id, case_id] ]
        germline_ranked_vcf // channel: [mandatory] [ val(group), val(meta), file(germline.ranked.vcf) ]

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

        // Fold the ranked germline calls (from GERMLINE_ANNOTATE) into the somatic case VCF,
        // joined by group, before the final filtering pass so it runs uniformly over both.
        // Profiles that disable germline reporting (params.germline == false, e.g. solid)
        // never populate germline_ranked_vcf, so skip straight to POST_ANNOTATION_FILTERS on
        // the somatic-only VCF instead of joining against an empty channel.
        if( params.germline ) {
            COMBINE_SOMATIC_GERMLINE (
                FIX_VEP.out.fixed_vcf
                    .join( germline_ranked_vcf )
                    .map { group, meta_somatic, somatic_vcf, meta_germline, germline_vcf ->
                        tuple(group, meta_somatic, somatic_vcf, germline_vcf)
                    }
            )
            ch_versions = ch_versions.mix(COMBINE_SOMATIC_GERMLINE.out.versions)
            ch_for_post_filter = COMBINE_SOMATIC_GERMLINE.out.vcf_combined
        }
        else {
            ch_for_post_filter = FIX_VEP.out.fixed_vcf
        }

        // add filters and mark duplicates
        POST_ANNOTATION_FILTERS { ch_for_post_filter }
        ch_versions = ch_versions.mix(POST_ANNOTATION_FILTERS.out.versions)

    emit:
        annotated_variants  =   FILTER_PANEL.out.vcf_panel                      // channel: [ val(group), val(vc), file(vcf.gz) ]
        finished_vcf        =   POST_ANNOTATION_FILTERS.out.filtered_vcf        // channel: [ val(group), val(vc), file(vcf.gz) ] (somatic+germline combined)
        versions            =   ch_versions                                     // channel: [ file(versions) ]
}
