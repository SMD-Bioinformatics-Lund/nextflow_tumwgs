include { GENS_VIZ              } from '../../modules/local/gens/main'
include { GENS_VIZ_V4              } from '../../modules/local/gens/main'

workflow VISUALIZE {
    take: 
        dnascope_vcf           
        gatk_count    

    main:
        ch_versions = Channel.empty()

        GENS_VIZ ( dnascope_vcf, gatk_count  )
        ch_versions = ch_versions.mix(GENS_VIZ.out.versions)

        GENS_VIZ_V4 ( GENS_VIZ.out.gens_for_v4  )


    emit:
        gens        =   GENS_VIZ.out.gens
        gens_v4     =   GENS_VIZ_V4.out.gens_v4
        versions    =   ch_versions                     // channel: [ file(versions) ]
}