process REHEADER_CRAM {
    label 'process_medium'
    label 'scratch'
    label 'stage'
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), file(cram), file(crai), file(bai)

    output:
		tuple val(group), val(meta), file("${out_cram_header}"), file("${out_cram_header}.crai"),   emit: cram_header_fixed
        path "versions.yml",                                                                        emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
	
        def args    = task.ext.args ?: ""
        out_cram_header = "${meta.id}.${meta.type}_chr.cram"

        """
        samtools view -H ${cram} \
        | sed 's/@SQ\\tSN:\\([0-9XYM]\\)/@SQ\\tSN:chr\\1/' \
        > header.chr.txt

        samtools reheader header.chr.txt ${cram} > ${out_cram_header}
        samtools index -@ ${task.cpus} ${out_cram_header}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(echo \$(samtools 2>&1) | sed 's/.*Version: //; s/ .*//')
        END_VERSIONS
        """
	
    stub:
        out_cram_header = "${meta.id}.${meta.type}_chr.cram"

        """
        touch ${out_cram_header} ${out_cram_header}.crai

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(echo \$(samtools 2>&1) | sed 's/.*Version: //; s/ .*//')
        END_VERSIONS
        """
}

process PELOPS_DUX4 {
    label 'process_medium'
    tag "${meta.id}"

    input:
		tuple val(group), val(meta), file(cram), file(crai)

    output:
    	tuple val(group), val(meta), path("*pelops.json"),  emit: pelops_dux4_json
        path "versions.yml",                                emit: versions
	
    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args ?: ""
        def prefix  = task.ext.prefix ?: "${meta.id}"    
        """
        pelops dux4r ${cram} \
            --json ${prefix}.pelops.json \
            --threads ${task.cpus} \
            ${args}
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            pelops: \$(echo \$(pelops version 2>&1))
        END_VERSIONS
        """

	stub:
        def args    = task.ext.args ?: ""
        def prefix  = task.ext.prefix ?: "${meta.id}" 
        """
        touch ${prefix}.pelops.json

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            pelops: \$(echo \$(pelops version 2>&1))
        END_VERSIONS
        """
}

process JSON_PELOPS_TO_VCF {
    label 'process_medium'
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), file(jsons)
    
    output:
    	tuple val(group), path ("*pelops.vcf"),         emit: pelops_dux4_vcf
        path "versions.yml",                            emit: versions
    
    when:
        task.ext.when == null || task.ext.when

    script:

        def args    = task.ext.args ?: ""
        def prefix  = task.ext.prefix ?: "${group}"

        if (meta.id.size() >= 2) {
            tumor_idx  = meta.type.findIndexOf { it == 'tumor' || it == 'T' }
            normal_idx = meta.type.findIndexOf { it == 'normal' || it == 'N' }
            id_tumor = meta[tumor_idx].id
            id_normal = meta[normal_idx].id
        
            """
            pelops_json_to_vcf.py --tumor-json ${jsons[tumor_idx]} --normal-json ${jsons[normal_idx]} --tumor-name ${id_tumor} --normal-name ${id_normal} --output ${prefix}.pelops.vcf
            
            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                python: \$(python --version 2>&1| sed -e 's/Python //g')
            END_VERSIONS
            """
        }
        else if (meta.id.size() == 1) {

            """
            pelops_json_to_vcf.py \\
                --tumor-json ${jsons[0]} \\ 
                --tumor-name ${meta[0].id} \\
                --output ${prefix}.pelops.vcf

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                python: \$(python --version 2>&1| sed -e 's/Python //g')
            END_VERSIONS
            """
        }

    stub:
        def args    = task.ext.args ?: ""
        def prefix  = task.ext.prefix ?: "${group}"

        if( meta.id.size() >= 2 ) {
            tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
            normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N' }
            id_tumor = meta[tumor_idx].id
            id_normal = meta[normal_idx].id 
            
            """
            echo $tumor_idx > log.txt
            echo $normal_idx >> log.txt
            echo $id_tumor  >> log.txt
            echo $id_normal >> log.txt

            touch ${prefix}.pelops.vcf
            cat <<-END_VERSIONS > versions.yml
            "${task.process}":    
                python: \$(python --version 2>&1| sed -e 's/Python //g')
            END_VERSIONS
            """
        }
        else if (meta.id.size() == 1) {
            """
            touch ${prefix}.pelops.vcf
            cat <<-END_VERSIONS > versions.yml
            "${task.process}":    
                python: \$(python --version 2>&1| sed -e 's/Python //g')
            END_VERSIONS
            """
        }
}








