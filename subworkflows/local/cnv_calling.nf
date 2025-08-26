#!/usr/bin/env nextflow

include { CNVKIT_BATCH                         } from '../../modules/local/cnvkit/main'
include { CNVKIT_GENS                          } from '../../modules/local/cnvkit/main'
include { CNVKIT_PLOT                          } from '../../modules/local/cnvkit/main'
include { CNVKIT_CALL                          } from '../../modules/local/cnvkit/main'
include { CNVKIT_CALL as CNVKIT_CALL_TC        } from '../../modules/local/cnvkit/main'
include { MERGE_GENS                           } from '../../modules/local/cnvkit/main'
include { CNVKIT_BATCH as CNVKIT_BACKBONE      } from '../../modules/local/cnvkit/main'
include { CNVKIT_BATCH as CNVKIT_EXONS         } from '../../modules/local/cnvkit/main'
include { GATKCOV_BAF                          } from '../../modules/local/GATK/main'
include { GATKCOV_COUNT                        } from '../../modules/local/GATK/main'
include { GATKCOV_CALL                         } from '../../modules/local/GATK/main'
include { GATK2VCF                             } from '../../modules/local/GATK/main'
include { GATK_COUNT_GERMLINE                  } from '../../modules/local/GATK/main'
include { GATK_CALL_PLOIDY                     } from '../../modules/local/GATK/main'
include { GATK_CALL_GERMLINE_CNV               } from '../../modules/local/GATK/main'
include { FILTER_MERGE_GATK                    } from '../../modules/local/GATK/main'
include { MERGE_GATK_TUMOR                     } from '../../modules/local/GATK/main'
include { POSTPROCESS                          } from '../../modules/local/GATK/main'
include { MANTA                                } from '../../modules/local/manta/main'
include { SVDB_MERGE_PANEL as JOIN_TUMOR       } from '../../modules/local/svdb/main'
include { SVDB_MERGE_PANEL as JOIN_NORMAL      } from '../../modules/local/svdb/main'
include { FILTER_MANTA as FILTER_MANTA_TUMOR   } from '../../modules/local/filters/main'
include { FILTER_MANTA as FILTER_MANTA_NORMAL  } from '../../modules/local/filters/main'


workflow CNV_CALLING {
    take: 
        bam_umi              // channel: [mandatory] [ val(group), val(meta), file(umi_bam), file(umi_bai), file(bqsr.table) ]
        germline_variants    // channel: [mandatory] [ val(group), file(vcf), file(tbi) ]
        meta                 // channel: [mandatory] [ [sample_id, group, sex, phenotype, paternal_id, maternal_id, case_id] ]
        bam_markdup          // channel: [mandatory] [ val(group), val(meta), file(dedup_bam), file(dedup_bai)]  ]
        gatk_ref             // channel: [mandatory] [ val(interger), val(part_of_genome) used for germline gatk-calling ]

    main:
        ch_versions = Channel.empty()
        ////////////////////////// CNVKIT /////////////////////////////////////////////////////
        // if backbone + exon pool differs in pool ratio do backbone and exons separatly
        if (!params.cnvkit_split) {
            CNVKIT_BATCH ( bam_umi, params.cnvkit_reference, "full" )
            batch_plot_cns = CNVKIT_BATCH.out.cnvkit_cns
            batch_plot_cnr = CNVKIT_BATCH.out.cnvkit_cnr
            ch_versions = ch_versions.mix(CNVKIT_BATCH.out.versions)

            // call, plot and export segments ::: cnvkit
            CNVKIT_PLOT ( batch_plot_cns.join(batch_plot_cnr, by:[0,1,3]).combine(germline_variants, by:[0]) )
            ch_versions = ch_versions.mix(CNVKIT_PLOT.out.versions)

            CNVKIT_CALL ( batch_plot_cns.join(batch_plot_cnr, by:[0,1,3]).combine(germline_variants, by:[0]), "false" )
            ch_versions = ch_versions.mix(CNVKIT_CALL.out.versions)

            CNVKIT_GENS ( batch_plot_cnr.combine(germline_variants, by:[0]) )
            ch_versions = ch_versions.mix(CNVKIT_GENS.out.versions)

            MERGE_GENS  ( CNVKIT_GENS.out.cnvkit_gens )
            ch_versions = ch_versions.mix(MERGE_GENS.out.versions)

            cnvkitplot = CNVKIT_PLOT.out.cnvkitplot
            cnvkit_hrd = CNVKIT_CALL.out.cnvkitsegment
            // tuple val(group), val(meta), val(part), file("${group}.${meta.id}.${meta.type}.${part}.vcf"), emit: cnvkit_vcf
            CNVKIT_VCF_TUMOR = CNVKIT_CALL.out.cnvkit_vcf.join(meta.filter( it -> it[1].type == "T" ) ).map{ val-> tuple(val[0], val[3], val[2] ) }
        }
        else {
            CNVKIT_BATCH ( bam_umi, params.cnvkit_reference, "full" )
            ch_versions = ch_versions.mix(CNVKIT_BATCH.out.versions)

            CNVKIT_EXONS ( bam_umi, params.cnvkit_reference_exons, "exons" )
            ch_versions = ch_versions.mix(CNVKIT_EXONS.out.versions)

            CNVKIT_BACKBONE ( bam_umi, params.cnvkit_reference_backbone, "backbone" )
            ch_versions = ch_versions.mix(CNVKIT_BACKBONE.out.versions)

            // call, plot and export segments ::: cnvkit
            CNVKIT_PLOT ( CNVKIT_BACKBONE.out.cnvkit_cns.join(CNVKIT_BACKBONE.out.cnvkit_cnr, by:[0,1,3]).combine(germline_variants, by:[0]) )
            ch_versions = ch_versions.mix(CNVKIT_PLOT.out.versions)

            // call without adjusting for purity //
            CNVKIT_CALL ( CNVKIT_EXONS.out.cnvkit_cns.join(CNVKIT_EXONS.out.cnvkit_cnr, by:[0,1,3])
                        .mix(
                            CNVKIT_BACKBONE.out.cnvkit_cns.join(CNVKIT_BACKBONE.out.cnvkit_cnr, by:[0,1,3]),
                            CNVKIT_BATCH.out.cnvkit_cns.join(CNVKIT_BATCH.out.cnvkit_cnr, by:[0,1,3])
                        )
                        .combine(germline_variants, by:[0]), "false" )
            ch_versions = ch_versions.mix(CNVKIT_CALL.out.versions)

            // call but adjust for purity, used for HRD //
            CNVKIT_CALL_TC( CNVKIT_BACKBONE.out.cnvkit_cns.join(CNVKIT_BACKBONE.out.cnvkit_cnr, by:[0,1,3]).combine(germline_variants, by:[0]), "true")
            ch_versions = ch_versions.mix(CNVKIT_CALL_TC.out.versions)

            // gens separatly for the two pools //
            CNVKIT_GENS ( CNVKIT_EXONS.out.cnvkit_cnr.mix(CNVKIT_BACKBONE.out.cnvkit_cnr).combine(germline_variants, by:[0]) )
            ch_versions = ch_versions.mix(CNVKIT_GENS.out.versions)

            // merge these //
            MERGE_GENS  ( CNVKIT_GENS.out.cnvkit_gens.groupTuple(by:[0,1]) )
            ch_versions = ch_versions.mix(MERGE_GENS.out.versions)

            // assign correct part full,exon,backbone to relevant upcoming analysis //
            cnvkitplot = CNVKIT_PLOT.out.cnvkitplot.filter { it -> it[2] == "backbone" }
            cnvkit_hrd = CNVKIT_CALL_TC.out.cnvkitsegment
            cnvkit_vcf = CNVKIT_CALL.out.cnvkit_vcf.filter { it -> it[1] == "full" }
            CNVKIT_VCF_TUMOR = cnvkit_vcf.join(meta.filter( it -> it[1].type == "T" ) ).map{ val-> tuple(val[0], val[3], val[2] ) }
        }

        ///////////////////////////////////////////////////////////////////////////////////////

        //////////////////////////// GATK SEGMENT CALLING /////////////////////////////////////
        // Do calling somatic CNV calling, use normal allelic counts for somatic as well     //
        ///////////////////////////////////////////////////////////////////////////////////////
        GATKCOV_BAF ( bam_umi )
        ch_versions = ch_versions.mix(GATKCOV_BAF.out.versions)

        GATKCOV_COUNT ( bam_umi )
        ch_versions = ch_versions.mix(GATKCOV_COUNT.out.versions)

        GATKCOV_CALL { GATKCOV_BAF.out.gatk_baf.join(GATKCOV_COUNT.out.gatk_count,by:[0,1]).groupTuple() }
        ch_versions = ch_versions.mix(GATKCOV_CALL.out.versions)

        GATK2VCF ( GATKCOV_CALL.out.gatcov_called.join(meta.filter( it -> it[1].type == "T" )) )
        ch_versions = ch_versions.mix(GATK2VCF.out.versions)

        MERGE_GATK_TUMOR ( GATK2VCF.out.tumor_vcf )
        ch_versions = ch_versions.mix(MERGE_GATK_TUMOR.out.versions)

        // Do germline calling for normal
        GATK_COUNT_GERMLINE ( bam_umi.filter { it -> it[1].type == "N" })
        ch_versions = ch_versions.mix(GATK_COUNT_GERMLINE.out.versions)

        GATK_CALL_PLOIDY ( GATK_COUNT_GERMLINE.out.count_germline )
        ch_versions = ch_versions.mix(GATK_CALL_PLOIDY.out.versions)

        GATK_CALL_GERMLINE_CNV( GATK_COUNT_GERMLINE.out.count_germline.join(GATK_CALL_PLOIDY.out.gatk_ploidy,by:[0,1]).groupTuple(by:[0,1]).combine(gatk_ref) )
        ch_versions = ch_versions.mix(GATK_CALL_GERMLINE_CNV.out.versions)

        CALLED = GATK_CALL_GERMLINE_CNV.out.gatk_call_germline.groupTuple(by:[0,1])
        PLOIDY = GATK_CALL_PLOIDY.out.gatk_ploidy.groupTuple(by:[0,1])

        POSTPROCESS ( CALLED.join(PLOIDY,by:[0,1]).combine(gatk_ref.groupTuple(by:[3])))
        ch_versions = ch_versions.mix(POSTPROCESS.out.versions)

        FILTER_MERGE_GATK ( POSTPROCESS.out.gatk_germline_segmentsvcf )
        ch_versions = ch_versions.mix(FILTER_MERGE_GATK.out.versions)

        /////////////////////////// MANTA /////////////////////////////////////////////////////
        MANTA ( bam_markdup.groupTuple(), params.bedgz, "CNV" )
        ch_versions = ch_versions.mix(MANTA.out.versions)

        // Join germline vcf
        // as manta is groupTupled, its output needs to be separated from meta, and rejoined or it wont be able to be mixable for SVDB
        MANTA_NORMAL = MANTA.out.manta_vcf_normal.join(meta.filter( it -> it[1].type == "N" ) ).map{ val-> tuple(val[0], val[2], val[1] ) }

        FILTER_MANTA_NORMAL(MANTA_NORMAL)
        ch_versions = ch_versions.mix(FILTER_MANTA_NORMAL.out.versions)

        GATK_NORMAL = FILTER_MERGE_GATK.out.gatk_normal_vcf.join(meta.filter( it -> it[1].type == "N" ) ).map{ val-> tuple(val[0], val[2], val[1] )}
        JOIN_NORMAL ( GATK_NORMAL.mix(FILTER_MANTA_NORMAL.out.filtered).groupTuple(by:[0,1]) )
        ch_versions = ch_versions.mix(JOIN_NORMAL.out.versions)

        // Join tumor vcf
        GATK_TUMOR = MERGE_GATK_TUMOR.out.tumor_vcf_merged
        // as manta is groupTupled, its output needs to be separated from meta, and rejoined or it wont be able to be mixable for SVDB
        MANTA_TUMOR = MANTA.out.manta_vcf_tumor.join(meta.filter( it -> it[1].type == "T" ) ).map{ val-> tuple(val[0], val[2], val[1] ) }

        FILTER_MANTA_TUMOR(MANTA_TUMOR)
        ch_versions = ch_versions.mix(FILTER_MANTA_TUMOR.out.versions)

        JOIN_TUMOR ( GATK_TUMOR.mix(FILTER_MANTA_TUMOR.out.filtered,CNVKIT_VCF_TUMOR).groupTuple(by:[0,1]) )
        ch_versions = ch_versions.mix(JOIN_TUMOR.out.versions)

    emit:
        gatcov_plot =   GATKCOV_CALL.out.gatcov_plot    // channel: [ val(group), file(modeled.png) ]
        cnvkit_plot =   cnvkitplot                      // channel: [ val(group), val(meta), val(part), file(cnvkit_overview.png) ]
        cnvkit_hrd  =   cnvkit_hrd                      // channel: [ val(group), val(meta), val(part), file(call.cns) ]
        tumor_vcf   =   JOIN_TUMOR.out.merged_vcf       // channel: [ val(group), val(vc), file(tumor.merged.vcf) ]
        normal_vcf  =   JOIN_NORMAL.out.merged_vcf      // channel: [ val(group), val(vc), file(normal.merged.vcf) ]
        gens        =   MERGE_GENS.out.dbload           // channel: [ val(group), val(meta), file(gens) ]
        versions    =   ch_versions                     // channel: [ file(versions) ]

}