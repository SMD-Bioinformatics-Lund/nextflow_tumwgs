process GENS_VIZ {
	label "process_low"
    tag "${meta.id}"
	
	input:
        tuple val(group), val(meta), path(vcf), path(tbi)
        tuple val(group), val(meta), path(cov_stand), path(cov_denoise)
		
	output:
		tuple val(group), val(meta), file("*baf.bed.gz*"), file("*cov.bed.gz*"), optional: true,    emit: gens
        tuple val(group), val(meta), file("*.gens"),                                                emit: dbload
        path "versions.yml",                                                                        emit: versions
    
    when:
        task.ext.when == null || task.ext.when
    
    script:
        process_group = group
        tumor_idx = 0
        if( meta.id.size() < 2 ) {
            process_group = group + '_TumorOnly'
            tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
        }
        def args    = task.ext.args   ?: ''

		"""
        generate_gens_data.pl ${cov_stand} ${vcf} ${meta.id} ${args}
        
        echo "gens load sample --sample-id ${meta.id} --case-id ${process_group} --genome-build 38 --baf ${params.gens_accessdir}/${meta.id}.baf.bed.gz --coverage ${params.gens_accessdir}/${meta.id}.cov.bed.gz" > ${meta.id}.gens
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
        END_VERSIONS
		"""

    stub:
        process_group = group
        tumor_idx = 0
        if( meta.id.size() < 2 ) {
            process_group = group + '_singel'
            tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
        }
        def args    = task.ext.args   ?: ''
        """
        touch ${meta.id}.cov.bed.gz 
        touch ${meta.id}.cov.bed.gz.tbi
        touch ${meta.id}.baf.bed.gz
        touch ${meta.id}.baf.bed.gz.tbi
        touch ${meta.id}.gens
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
        END_VERSIONS
        """
}