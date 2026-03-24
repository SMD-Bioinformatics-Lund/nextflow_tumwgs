#!/usr/bin/env nextflow

include { REHEADER_CRAM                          } from '../../modules/local/dux4_igh/main'
include { PELOPS_DUX4                            } from '../../modules/local/dux4_igh/main'
include { JSON_PELOPS_TO_VCF                     } from '../../modules/local/dux4_igh/main'

workflow DUX4IGH_CALLING {
    take: 
        cram_dedup            // channel: [mandatory] [ val(group), val(meta), file(cram), file(crai), file(bai) ]
        meta                 // channel: [mandatory] [ [sample_id, group, sex, phenotype, paternal_id, maternal_id, case_id] ]

    main:
        ch_versions = Channel.empty()

        //////////////////////////// IGH - DUX4 calling/////////////////////////////////////
        // DUX4-IGH fusion calling is based on Pelops, which requires the CRAM header to have "chr" prefix for chromosome names.
        ///////////////////////////////////////////////////////////////////////////////////////
        REHEADER_CRAM(  cram_dedup  )
        ch_versions = ch_versions.mix(REHEADER_CRAM.out.versions)

        PELOPS_DUX4 ( REHEADER_CRAM.out.cram_header_fixed)
        ch_versions = ch_versions.mix(PELOPS_DUX4.out.versions)
        
        PELOPS_DUX4.out.pelops_dux4_json.groupTuple().view { println "PELOPS_DUX4 output: ${it}" }

        JSON_PELOPS_TO_VCF ( PELOPS_DUX4.out.pelops_dux4_json.groupTuple() )
        ch_versions = ch_versions.mix(JSON_PELOPS_TO_VCF.out.versions)  

    emit:
        pelops_dux4_vcf   = JSON_PELOPS_TO_VCF.out.pelops_dux4_vcf
        versions    =   ch_versions     // channel: [ file(versions) ]

}