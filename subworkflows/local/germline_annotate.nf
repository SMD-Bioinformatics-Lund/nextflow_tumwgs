#!/usr/bin/env nextflow

include { NORMALIZE_VCF            } from '../../modules/local/filters/main'
include { INTERSECT_CODING         } from '../../modules/local/filters/main'
include { ANNOTATE_VEP             } from '../../modules/local/filters/main'
include { GERMLINE_EVALUATE        } from '../../modules/local/germline_evaluate/main'

workflow GERMLINE_ANNOTATE {
    take:
        germline_combined_vcf   // channel: [mandatory] [ val(group), val(meta), file(vcf.gz), file(vcf.gz.tbi) ]

    main:
        ch_versions = Channel.empty()

        NORMALIZE_VCF { germline_combined_vcf }
        ch_versions = ch_versions.mix(NORMALIZE_VCF.out.versions)

        INTERSECT_CODING ( NORMALIZE_VCF.out.vcf_norm, params.gene_regions )
        ch_versions = ch_versions.mix(INTERSECT_CODING.out.versions)

        ANNOTATE_VEP { INTERSECT_CODING.out.vcf_intersected }
        ch_versions = ch_versions.mix(ANNOTATE_VEP.out.versions)

        GERMLINE_EVALUATE { ANNOTATE_VEP.out.vcf_vep }
        ch_versions = ch_versions.mix(GERMLINE_EVALUATE.out.versions)

    emit:
        ranked_vcf = GERMLINE_EVALUATE.out.vcf_ranked      // channel: [ val(group), val(meta), file(vcf) ]
        versions   = ch_versions                            // channel: [ file(versions) ]
}
