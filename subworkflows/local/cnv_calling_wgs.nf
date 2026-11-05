#!/usr/bin/env nextflow

include { GATKCOV_BAF                          } from '../../modules/local/GATK/main'
include { GATKCOV_COUNT                        } from '../../modules/local/GATK/main'
include { GATKCOV_CALL                         } from '../../modules/local/GATK/main'
include { GATKCOV_CALL_GERMLINE                } from '../../modules/local/GATK/main'
include { OVERLAP_GENES                        } from '../../modules/local/filters/main'
include { FILTER_CNVS_PANEL                    } from '../../modules/local/filters/main'
include { ASCAT_PREPARE                        } from '../../modules/local/ascat/main'
include { ASCAT_FIT                            } from '../../modules/local/ascat/main'
include { ASCAT_LOH_HEH                        } from '../../modules/local/ascat/main'


workflow CNV_CALLING {
    take:
        cram_dedup            // channel: [mandatory] [ val(group), val(meta), file(cram), file(crai), file(bai) ]
        meta                 // channel: [mandatory] [ [sample_id, group, sex, phenotype, paternal_id, maternal_id, case_id] ]

    main:
        ch_versions = Channel.empty()

        //////////////////////////// GATK SEGMENT CALLING /////////////////////////////////////
        // Do calling somatic CNV calling, use normal allelic counts for somatic as well     //
        ///////////////////////////////////////////////////////////////////////////////////////
        GATKCOV_BAF (  cram_dedup  )
        ch_versions = ch_versions.mix(GATKCOV_BAF.out.versions)

        GATKCOV_COUNT ( cram_dedup )
        ch_versions = ch_versions.mix(GATKCOV_COUNT.out.versions)

        ch_gatkcov_grouped = GATKCOV_BAF.out.gatk_baf.join(GATKCOV_COUNT.out.gatk_count,by:[0,1]).groupTuple()

        GATKCOV_CALL { ch_gatkcov_grouped }
        ch_versions = ch_versions.mix(GATKCOV_CALL.out.versions)

        // GATKCOV_CALL_GERMLINE calls segments on the normal sample; skip groups
        // that have no normal (tumor-only) since there is nothing to call.
        GATKCOV_CALL_GERMLINE { ch_gatkcov_grouped.filter { grp, sample_meta, allele, stdCR, denoised -> sample_meta.size() == 2 } }
        ch_versions = ch_versions.mix(GATKCOV_CALL.out.versions)

        OVERLAP_GENES { GATKCOV_CALL.out.gatkcov_called }
        ch_versions = ch_versions.mix(OVERLAP_GENES.out.versions)

        FILTER_CNVS_PANEL { OVERLAP_GENES.out.annotated_bed }
        ch_versions = ch_versions.mix(OVERLAP_GENES.out.versions)

        ////////////////////////////// ASCAT LOH / HeH ////////////////////////////////////////
        // Allele-specific copy number (nMajor/nMinor) for tumor/normal pairs only. ASCAT    //
        // needs the sex; groups whose normal has no F/M sex in the samplesheet are skipped.  //
        ///////////////////////////////////////////////////////////////////////////////////////
        ch_ascat_in = cram_dedup
            .groupTuple()
            .filter { grp, sample_meta, cram, crai, bai -> sample_meta.size() == 2 }
            .filter { grp, sample_meta, cram, crai, bai ->
                def n = sample_meta.findIndexOf { it.type == 'normal' || it.type == 'N' }
                def ok = n >= 0 && "${sample_meta[n].sex}".toUpperCase() in ['F', 'M', 'FEMALE', 'MALE', 'XX', 'XY']
                if (!ok) log.warn "[ASCAT] ${grp}: normal sample has no F/M sex — skipping ASCAT"
                ok
            }

        ASCAT_PREPARE { ch_ascat_in }
        ch_versions = ch_versions.mix(ASCAT_PREPARE.out.versions)

        ASCAT_FIT { ASCAT_PREPARE.out.logr_baf }
        ch_versions = ch_versions.mix(ASCAT_FIT.out.versions)

        ASCAT_LOH_HEH { ASCAT_FIT.out.segments }
        ch_versions = ch_versions.mix(ASCAT_LOH_HEH.out.versions)

    emit:
        tum_plot    =   GATKCOV_CALL.out.gatkcov_plot
        norm_plot   =   GATKCOV_CALL_GERMLINE.out.gatkcov_germlline_plot
        count       =   GATKCOV_COUNT.out.gatk_count
        bed         =   FILTER_CNVS_PANEL.out.vcf_panel_bed
        json        =   FILTER_CNVS_PANEL.out.vcf_panel_json
        ascat_segments = ASCAT_FIT.out.segments             // [ group, meta, purity_ploidy.tsv, segments.tsv, segments_raw.tsv ]
        ascat_loh_heh  = ASCAT_LOH_HEH.out.loh_heh          // [ group, meta, loh.tsv, subclonal.tsv, gains.tsv, chrom.tsv, heh.txt ]
        ascat_cn_summary = ASCAT_LOH_HEH.out.cn_summary     // [ group, meta, cn_summary.json, arms.tsv, events.tsv ] (standardised, for reporting)
        versions    =   ch_versions                     // channel: [ file(versions) ]

}
