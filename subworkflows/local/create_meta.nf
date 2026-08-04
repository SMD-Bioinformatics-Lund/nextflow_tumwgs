#!/usr/bin/env nextflow

// might need to add a check to csv? //
include { CSV_CHECK      } from '../../modules/local/check_input/main'

workflow CHECK_INPUT {
    take:
        csv     // file(csv)
        paired  // boolean

    main:
        CSV_CHECK(csv)
        checkedCSV = CSV_CHECK.out.csv.splitCsv(header:true, sep:',').set { csvmap }

		group_size = csvmap.map { [it.group, it] }.groupTuple().map { group, rows -> [group, rows.size()] }

		fastq     = csvmap.map { create_fastq_channel(it) }
		meta      = csvmap.map { [it.group, it] }
						.combine(group_size, by: 0)
						.map { group, row, size -> create_samples_channel(row, size) }

    emit:
        fastq
        bam
        cram
        vcf
        meta
}


// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channel(LinkedHashMap row, paired) {
	// create meta map
	def meta = [:]
	meta.id             	= row.id
	meta.group              = row.group
	meta.diagnosis          = row.diagnosis
	meta.type               = row.type
    meta.platform			= row.platform
	meta.clarity_sample_id  = row.clarity_sample_id
	meta.ffpe               = row.containsKey("ffpe") && row.ffpe ? true : false
	meta.purity             = (row.containsKey("purity") ? row.purity : false)
	meta.sequencing_run     = row.sequencing_run
	meta.reads              = (row.containsKey("n_reads") ? row.n_reads : false)
    meta.sex                = (row.containsKey("sex") ? row.sex : false)
	meta.clarity_pool_id    = row.clarity_pool_id
    meta.paired             = paired
    def sub = false
    if (meta.reads && params.sample) {  
        if (meta.reads.toInteger() > params.sample_val) {
            sub = (params.sample_val / meta.reads.toInteger()).round(2)
            if (sub == 1.00) sub = 0.99
        }
    }
    meta.sub = sub
	// add path(s) of the fastq file(s) to the meta map
	def fastq_meta = []
	fastq_meta = [row.group, meta, file(row.read1), file(row.read2) ]

	return fastq_meta
}

// Function to get a list of metadata (e.g. pedigree, case id) from the sample; [ meta ]
def create_samples_channel(LinkedHashMap row, int group_size) {
	def meta                = [:]
	meta.id                 = row.id
	meta.group				= row.group
	meta.diagnosis          = row.diagnosis
	meta.sex				= row.sex
	meta.type               = row.type
	meta.tumor_only         = (group_size == 1)
	meta.clarity_sample_id  = row.clarity_sample_id
	meta.ffpe               = (row.containsKey("ffpe") ? row.ffpe : false)
	meta.purity             = (row.containsKey("purity") ? row.purity : false)
	meta.sequencing_run     = row.sequencing_run
	meta.reads              = (row.containsKey("n_reads") ? row.n_reads : false)
	meta.clarity_pool_id    = row.clarity_pool_id
	sub = false
	if (meta.reads && params.sample) {  
		if (meta.reads.toInteger() > params.sample_val) {
			sub = (params.sample_val / meta.reads.toInteger()).round(2)
			if (sub == 1.00){
				sub = 0.99
			}
		}
		else {
			sub = false
		}
	}
	meta.sub = sub

	def sample_meta = [row.group, meta]
	return sample_meta
}
