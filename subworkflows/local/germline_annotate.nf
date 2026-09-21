#!/usr/bin/env nextflow

include { NORMALIZE_VCF            } from '../../modules/local/filters/main'
include { INTERSECT_CODING         } from '../../modules/local/filters/main'
include { FORMAT_GERMLINE_VCF      } from '../../modules/local/filters/main'
include { ANNOTATE_VEP             } from '../../modules/local/filters/main'
include { MARK_GERMLINES           } from '../../modules/local/filters/main'
include { SELECT_GERMLINE          } from '../../modules/local/filters/main'

workflow GERMLINE_ANNOTATE {
    take:
        germline_combined_vcf   // channel: [mandatory] [ val(group), val(meta), file(vcf.gz), file(vcf.gz.tbi) ]

    main:
        ch_versions = Channel.empty()

        NORMALIZE_VCF { germline_combined_vcf }
        ch_versions = ch_versions.mix(NORMALIZE_VCF.out.versions)

        INTERSECT_CODING ( NORMALIZE_VCF.out.vcf_norm, params.gene_regions )
        ch_versions = ch_versions.mix(INTERSECT_CODING.out.versions)

        // DNAscope GT:AD:DP:GQ:PL -> GT:VAF:VD:DP, the format mark_germlines.pl expects
        FORMAT_GERMLINE_VCF { INTERSECT_CODING.out.vcf_intersected }
        ch_versions = ch_versions.mix(FORMAT_GERMLINE_VCF.out.versions)

        ANNOTATE_VEP { FORMAT_GERMLINE_VCF.out.vcf_agg }
        ch_versions = ch_versions.mix(ANNOTATE_VEP.out.versions)

        // flags variants that are in the assay genes and present in the normal as FILTER=GERMLINE
        MARK_GERMLINES { ANNOTATE_VEP.out.vcf_vep }
        ch_versions = ch_versions.mix(MARK_GERMLINES.out.versions)

        // keep only the flagged records, they are what gets merged into the somatic VCF
        SELECT_GERMLINE { MARK_GERMLINES.out.vcf_germline }
        ch_versions = ch_versions.mix(SELECT_GERMLINE.out.versions)

    emit:
        ranked_vcf = SELECT_GERMLINE.out.vcf_germline      // channel: [ val(group), val(meta), file(vcf) ]
        versions   = ch_versions                            // channel: [ file(versions) ]
}
