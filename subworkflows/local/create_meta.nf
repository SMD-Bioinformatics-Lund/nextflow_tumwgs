#!/usr/bin/env nextflow

// might need to add a check to csv? //
include { CSV_CHECK      } from '../../modules/local/check_input/main'

workflow CHECK_INPUT {
	take:
		csv		// file(csv)

	main:
		CSV_CHECK ( csv )
		checkedCsv = CSV_CHECK.out.csv.splitCsv( header:true, sep:',').set { csvmap }

		fastq     = csvmap.map { create_fastq_channel(it) }
		meta      = csvmap.map { create_samples_channel(it) }

	emit:
		fastq        // channel: [ val(meta), [ reads ] ]
		meta         // channel: [ sample_id, sex, phenotype, paternal_id, maternal_id, case_id ]

}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channel(LinkedHashMap row) {
	// create meta map
	def meta = [:]
	meta.id					= row.id
	meta.group              = row.group
	meta.diagnosis          = row.diagnosis
	meta.sex				= row.sex
	meta.type               = row.type
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
	// add path(s) of the fastq file(s) to the meta map
	def fastq_meta = []
	fastq_meta = [row.group, meta, file(row.read1), file(row.read2) ]

	return fastq_meta
}

// Function to get a list of metadata (e.g. pedigree, case id) from the sample; [ meta ]
def create_samples_channel(LinkedHashMap row) {
	def meta                = [:]
	meta.id                 = row.id
	meta.group				= row.group
	meta.diagnosis          = row.diagnosis
	meta.sex				= row.sex
	meta.type               = row.type
	meta.clarity_sample_id  = row.clarity_sample_id
	meta.ffpe               = (row.containsKey("ffpe") ? row.ffpe : false)
	meta.purity             = (row.containsKey("purity") ? row.purity : false)
	meta.sequencing_run     = row.sequencing_run
	meta.reads              = (row.containsKey("n_reads") ? row.n_reads : false)
	meta.clarity_pool_id    = row.clarity_pool_id

	def sample_meta = [row.group, meta]
	return sample_meta
}
