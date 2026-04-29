#!/usr/bin/env nextflow

include { SENTIEON_QC          } from '../../modules/local/sentieon/main'
include { COLLECT_QC           } from '../../modules/local/sentieon/main'
include { QC_TO_CDM            } from '../../modules/local/qc/main'
include { ALLELE_CALL          } from '../../modules/local/idSnp/main'
include { SNP_CHECK            } from '../../modules/local/idSnp/main'
include { PAIRGEN_CDM          } from '../../modules/local/idSnp/main'
include { VERIFYBAMID2         } from '../../modules/local/verfifybam2/main'
include { MERGE_QC_JSON        } from '../../modules/local/qc/main'
include { CREATE_PED_FILES     } from '../../modules/local/somalier/main'
include { SOMALIER_QC          } from '../../modules/local/somalier/main'



workflow BAM_QC {
    take:
        bam_bqsr        // channel: [ val(group), val(meta), file(bam), file(bai) ]
        bam_dedup       // channel: [ val(group), val(meta), file(cram), file(crai), file(bai) ]
        dedup_metrics   // channel: [ val(group), val(meta), file(dedup_metrics) ]
        meta       // channel: [ val(group), val(meta)]

    main:
        ch_versions = Channel.empty()
        ch_qc_json  = Channel.empty()
        
        SENTIEON_QC ( bam_dedup.join(dedup_metrics, by:[0,1]) )
        ch_versions = ch_versions.mix(SENTIEON_QC.out.versions)

        COLLECT_QC ( SENTIEON_QC.out.qc )
        ch_qc_json = ch_qc_json.mix(COLLECT_QC.out.alignment_qc)
        ch_versions = ch_versions.mix(COLLECT_QC.out.versions)

        VERIFYBAMID2 (bam_dedup)
        ch_versions = ch_versions.mix(VERIFYBAMID2.out.versions)
        ch_qc_json  = ch_qc_json.join(VERIFYBAMID2.out.contamination_json, by:[0,1])

        //ch_qc_json.view()

        CREATE_PED_FILES ( meta )
        ch_versions = ch_versions.mix(CREATE_PED_FILES.out.versions)
        
        ch_PED = CREATE_PED_FILES.out.ped_file 

        ch_PED_grouped = ch_PED
            .groupTuple(by: 0)
            .map { group_id, meta_list, ped_list ->
                def meta_by_type = [:]
                def ped_by_type = [:]
        
                meta_list.eachWithIndex { meta, idx ->
                meta_by_type[meta.type] = meta
                ped_by_type[meta.type] = ped_list[idx]
            }
        
            // Order N then T
            def ordered_types = ['N', 'T']
            def ordered_meta = ordered_types.collect { meta_by_type[it] }
            def ordered_peds = ordered_types.collect { ped_by_type[it] }
        
            [group_id, ordered_meta, ordered_peds]
        }
        ch_BAM_PED_grouped = bam_dedup.groupTuple().join(ch_PED_grouped, by: 0)

        ch_somalier = ch_BAM_PED_grouped.map { group_id, bam_meta, crams, crais, bais, ped_meta, peds ->
            [group_id, bam_meta, crams, crais, bais, peds]
        }
        
        ch_somalier.view { println "Somalier input: ${it}"  }
        SOMALIER_QC ( ch_somalier )
        ch_versions = ch_versions.mix(SOMALIER_QC.out.versions)

        MERGE_QC_JSON (ch_qc_json)
        ch_versions = ch_versions.mix(MERGE_QC_JSON.out.versions)

        QC_TO_CDM ( MERGE_QC_JSON.out.qc_json )

        // Check genotypes of ID-SNPs
        ALLELE_CALL (bam_bqsr)
        ch_versions = ch_versions.mix(ALLELE_CALL.out.versions)

        SNP_CHECK(ALLELE_CALL.out.sample_id_genotypes.groupTuple())
        ch_versions = ch_versions.mix(SNP_CHECK.out.versions)

        PAIRGEN_CDM (SNP_CHECK.out.idsnp_checked)

    emit:

        qcdone                  =   QC_TO_CDM.out.cdm_done                        // channel: [ val(group), val(meta), file
        versions                =   ch_versions                                  // channel: [ file(versions) ]
        dedup_cram_is_metrics   =   SENTIEON_QC.out.dedup_cram_is_metrices      // [ val(group), val(meta), file(dedup_cram), file(dedup_crai), file(dedup_bai) file(is_metrics) ]
        //qcdone                  =   COLLECT_QC.out.qc_cdm   
}