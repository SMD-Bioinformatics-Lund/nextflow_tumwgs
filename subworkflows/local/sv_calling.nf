#!/usr/bin/env nextflow

include { MANTA                                } from '../../modules/local/manta/main'
include { MANTA_SV                             } from '../../modules/local/manta/main'
include { SNPEFF                               } from '../../modules/local/snpeff/main'
include { SNPEFF_SV_ANN                        } from '../../modules/local/snpeff/main'
include { COMBINE_FUSIONS                      } from '../../modules/local/snpeff/main'
include { FILTER_FUSIONS_PANEL                 } from '../../modules/local/filters/main'


workflow SV_CALLING {
    take: 
        cram_dedup            // channel: [mandatory] [ val(group), val(meta), file(cram), file(crai), file(bai) ]
        meta                 // channel: [mandatory] [ [sample_id, group, sex, phenotype, paternal_id, maternal_id, case_id] ]

    main:
        ch_versions = Channel.empty()

        //////////////////////////// GATK SEGMENT CALLING /////////////////////////////////////
        // Do calling somatic CNV calling, use normal allelic counts for somatic as well     //
        ///////////////////////////////////////////////////////////////////////////////////////

        MANTA { cram_dedup }
        ch_versions = ch_versions.mix(MANTA.out.versions)

        SNPEFF { MANTA.out.manta_vcf_tumor }
        ch_versions = ch_versions.mix(SNPEFF.out.versions)

        MANTA_SV { MANTA.out.manta_vcf_tumor }
        ch_versions = ch_versions.mix(MANTA_SV.out.versions)
        
        SNPEFF_SV_ANN { MANTA_SV.out.inv_vcf_tumor }
        ch_versions = ch_versions.mix(SNPEFF_SV_ANN.out.versions)

        FILTER_FUSIONS_PANEL { SNPEFF.out.snpeff_vcf }
        ch_versions = ch_versions.mix( FILTER_FUSIONS_PANEL.out.versions)

        //ch_cmdfusion = FILTER_FUSIONS_PANEL.out.sv_panel.combine(SNPEFF_SV_ANN.out.snpeff_CMD)
        //ch_cmdfusion.view()
        //COMBINE_FUSIONS { ch_cmdfusion }

        COMBINE_FUSIONS  ( FILTER_FUSIONS_PANEL.out.sv_panel, SNPEFF_SV_ANN.out.snpeff_BND_TANDEM )
        ch_versions = ch_versions.mix( COMBINE_FUSIONS.out.versions)

    emit:
        // fusions     =   FILTER_FUSIONS_PANEL.out.sv_panel
        fusions     =   COMBINE_FUSIONS.out.sv_CMD
        versions    =   ch_versions                     // channel: [ file(versions) ]

}