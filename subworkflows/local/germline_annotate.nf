#!/usr/bin/env nextflow

include { ANNOTATE_VEP             } from '../../modules/local/filters/main'
include { GERMLINE_EVALUATE        } from '../../modules/local/germline_evaluate/main'

workflow GERMLINE_ANNOTATE {
    take:
        germline_combined_vcf   // channel: [mandatory] [ val(group), val(meta), file(vcf.gz), file(vcf.gz.tbi) ]

    main:
        ch_versions = Channel.empty()

        ANNOTATE_VEP { germline_combined_vcf.map{ group, meta, vcf, tbi -> tuple(group, meta, vcf) } }
        ch_versions = ch_versions.mix(ANNOTATE_VEP.out.versions)

        GERMLINE_EVALUATE { ANNOTATE_VEP.out.vcf_vep }
        ch_versions = ch_versions.mix(GERMLINE_EVALUATE.out.versions)

    emit:
        ranked_vcf = GERMLINE_EVALUATE.out.vcf_ranked      // channel: [ val(group), val(meta), file(vcf) ]
        versions   = ch_versions                            // channel: [ file(versions) ]
}
