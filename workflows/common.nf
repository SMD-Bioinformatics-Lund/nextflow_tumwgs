nextflow.enable.dsl = 2

include { CHECK_INPUT                   } from '../subworkflows/local/create_meta'
include { SAMPLE                        } from '../subworkflows/local/sample'
include { CUSTOM_DUMPSOFTWAREVERSIONS   } from '../modules/local/custom/dumpsoftwareversions/main'
include { ALIGN_SENTIEON                } from '../subworkflows/local/align_sentieon'
include { BAM_QC                        } from '../subworkflows/local/bam_qc'
include { SNV_CALLING                   } from '../subworkflows/local/snv_calling'
include { SNV_ANNOTATE                  } from '../subworkflows/local/snv_annotate'
include { CNV_CALLING                   } from '../subworkflows/local/cnv_calling_wgs'
include { SV_CALLING                    } from '../subworkflows/local/sv_calling'
include { VISUALIZE                     } from '../subworkflows/local/visualize'
include { ADD_TO_DB                     } from '../subworkflows/local/add_to_db'


csv = file(params.csv)

Channel
    .from( 0..params.bwa_shards-1)
    .set { ch_bwa_shards }

Channel
    .fromPath(params.genomic_shards_file)
    .splitCsv(header:false)
    .set {ch_shards }

Channel
    .fromPath("${params.regions_bed}")
    .ifEmpty { exit 1, "Regions bed file not found: ${params.regions_bed}" }
    .splitText( by: 20000, file: 'bedpart.bed' )
    .set { ch_beds }

workflow SWGP_COMMON {
    
    ch_versions = Channel.empty()

    // Checks input, creates meta-channel and decides whether data should be downsampled //
    CHECK_INPUT ( Channel.fromPath(csv) )

    // Downsample if meta.sub == value and not false //
    SAMPLE ( CHECK_INPUT.out.fastq )  
    .set{ ch_trim }
    ch_versions = ch_versions.mix(ch_trim.versions)

    // Do alignment if downsample was false and mix with SAMPLE subworkflow output
    ALIGN_SENTIEON ( 
        ch_bwa_shards,
        ch_trim.fastq_trim,
        CHECK_INPUT.out.meta
    )
    .set { ch_mapped }
    ch_versions = ch_versions.mix(ch_mapped.versions)

    BAM_QC (
        ch_mapped.bam_bqsr,
        ch_mapped.cram_dedup,
        ch_mapped.dedup_metrics
    )
    .set { ch_qc }
    ch_versions = ch_versions.mix(ch_qc.versions)

    SNV_CALLING ( 
        ch_mapped.bam_bqsr.groupTuple(),
        ch_mapped.cram_dedup,
        ch_beds,
        CHECK_INPUT.out.meta,
    )
    .set { ch_vcf }
    ch_versions = ch_versions.mix(ch_vcf.versions)


    SNV_ANNOTATE (
        ch_vcf.agg_vcf,
        ch_vcf.concat_vcfs,
        CHECK_INPUT.out.meta
    )
    .set { ch_vcf_anno }
    ch_versions = ch_versions.mix(ch_vcf_anno.versions)


    CNV_CALLING ( 
        ch_mapped.cram_dedup, 
        CHECK_INPUT.out.meta
    )
    .set { ch_cnvcalled }
    ch_versions = ch_versions.mix(ch_cnvcalled.versions)
    
    
    SV_CALLING (
                ch_mapped.cram_dedup.groupTuple(), 
                CHECK_INPUT.out.meta,
    )
    .set { ch_svcalled }
    ch_versions = ch_versions.mix(ch_svcalled.versions)
    
    VISUALIZE (
                ch_vcf.dnascope_vcf, 
                ch_cnvcalled.count,
    )
    .set { ch_visualize }
    ch_versions = ch_versions.mix(ch_visualize.versions)
    
    ADD_TO_DB (
        ch_vcf_anno.finished_vcf,
        ch_cnvcalled.vcf,
        ch_svcalled.fusions,
        ch_cnvcalled.tum_plot,
    )
    
    
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml'),
        CHECK_INPUT.out.meta
    )
}

workflow.onComplete {

    def msg = """\
        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        scriptFile  : ${workflow.scriptFile}
        workDir     : ${workflow.workDir}
        csv         : ${params.csv}
        exit status : ${workflow.exitStatus}
        errorMessage: ${workflow.errorMessage}
        errorReport :
        """
        .stripIndent()

    def error = """\
        ${workflow.errorReport}
        """
        .stripIndent()

    base = csv.getBaseName()
    logFile = file("${params.resultsdir}/cron/logs/" + base + ".complete")
    logFile.text = msg
    logFile.append(error)
}
