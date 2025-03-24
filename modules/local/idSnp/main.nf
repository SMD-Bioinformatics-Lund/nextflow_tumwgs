process ALLELE_CALL {
    label 'process_single'
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), file(bam), file(bai)

    output:
        tuple val(group), val(meta), file("*final.vcf"), file("*genotypes.json"),       emit:   sample_id_genotypes
        path "versions.yml",                                                            emit:   versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix  = task.ext.prefix ?: "${meta.id}"
        def args    = task.ext.args  ?: ""
        def args2   = task.ext.args2 ?: ""
        def args3   = task.ext.args3 ?: ""
        """
        bcftools mpileup $args $bam | bcftools call $args2   > ${prefix}.raw.vcf
        bcftools annotate $args3 -o ${prefix}.final.vcf ${prefix}.raw.vcf
        bcftools query -f '%ID\\t[%GT]\\n' ${prefix}.final.vcf > ${prefix}.genotypes
        genotype2json.py ${prefix}.genotypes ${prefix}.genotypes.json

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bcftools: \$(echo \$(bcftools --version 2>&1) | sed 's/bcftools //; s/ .*//')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        touch ${prefix}.final.vcf
        touch ${prefix}.genotypes.json

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bcftools: \$(echo \$(bcftools --version 2>&1) | sed 's/bcftools //; s/ .*//')
        END_VERSIONS
        """
}

process SNP_CHECK {
    label 'process_single'
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), file(vcfs), file(genotypes)

    output:
        tuple val(group), val(meta), file("*.idsnp"),           emit: idsnp_checked
        path "versions.yml", optional: true,                    emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        
        def args                    = task.ext.args  ?: ""
        tumor_idx                   = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
        normal_idx                  = meta.type.findIndexOf{ it == 'normal' || it == 'N' }
        normal_id                   = meta.id[normal_idx]
        tumor_id                    = meta.id[tumor_idx]
                    
        tumor_seq_run               = meta.sequencing_run[tumor_idx]
        normal_seq_run               = meta.sequencing_run[normal_idx]

	    normalvcf                   = vcfs[normal_idx]
	    tumorvcf                    = vcfs[tumor_idx]

        norGTjson                   = genotypes[normal_idx]
        tumGTjson                   = genotypes[tumor_idx]

        if(meta.id.size() == 2) {
            """
            idsnp_controller.pl \\
                --vcf_sample $tumorvcf  \\
                --vcf_control $normalvcf \\
                --sample  $tumor_id \\
                --control $normal_id \\
                $args

            cp s${tumor_id}_c${normal_id}.json  ${tumor_id}.json
            cp s${tumor_id}_c${normal_id}.json  ${normal_id}.json
            rm s${tumor_id}_c${normal_id}.json
            today_date=\$(date -u +"%Y-%m-%dT%H:%M:%SZ")

            echo '{"partner" : "${tumor_id}","sequencing_run" : "${tumor_seq_run}","analysis_date" : "'\${today_date}'" }' > ${normal_id}_partner_info.json
            echo '{"partner" : "${normal_id}","sequencing_run" : "${normal_seq_run}","analysis_date" : "'\${today_date}'" }' > ${tumor_id}_partner_info.json
            
            combinejsons.py  ${tumor_id}.json ${tumGTjson} --partner_run_json_file  ${tumor_id}_partner_info.json ${tumor_id}.idsnp
            combinejsons.py  ${normal_id}.json ${norGTjson} --partner_run_json_file ${normal_id}_partner_info.json ${normal_id}.idsnp
        
            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
            END_VERSIONS
            """
        } else {
            """
            echo "Not applicable" > s${tumor_id}.csv
            echo  '{ "is_paired_sample" : false}' > ${tumor_id}.json
            combinejsons.py  ${tumor_id}.json ${tumGTjson} ${tumor_id}.idsnp
            
            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
            END_VERSIONS
            """
        }

    stub:
        tumor_idx   = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
        normal_idx  = meta.type.findIndexOf{ it == 'normal' || it == 'N' }
        normal_id   = meta.id[normal_idx]
        tumor_id    = meta.id[tumor_idx]
        
        if(meta.id.size() == 2) {
            """ 
            touch s${tumor_id}_c${normal_id}.csv
            touch ${tumor_id}.idsnp
            touch ${normal_id}.idsnp


            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
            END_VERSIONS
            """
        } else {
        
            """
            echo "Not applicable" > s${tumor_id}.csv
            echo  '{ "is_paired_sample" : false }' > ${tumor_id}.json
            touch ${tumor_id}.idsnp

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
            END_VERSIONS
            """
        }
}

process PAIRGEN_CDM {
    label 'process_single'
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), file(jsons)

    output:
        tuple val(group), val(meta), file("*.pairgen"), emit: isnsp_cdm_done

    when:
        task.ext.when == null || task.ext.when

    script:
      def args    = task.ext.args  ?: ""
      if(meta.id.size() == 2) {
        
        tumor_idx   = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
        normal_idx  = meta.type.findIndexOf{ it == 'normal' || it == 'N' }
        normal_id   = meta.id[normal_idx]
        tumor_id    = meta.id[tumor_idx]
	    normaljson  = jsons[normal_idx]
	    tumorjson    = jsons[tumor_idx]

        """
        echo "--overwrite --sample-id ${meta.id[tumor_idx]} --sequencing-run ${meta.sequencing_run[tumor_idx]} --assay ${params.cdm} --id-snp ${params.outdir}/${params.subdir}/QC/${tumorjson} " > ${meta.id[tumor_idx]}.pairgen

             
        echo "--overwrite --sample-id ${meta.id[normal_idx]} --sequencing-run ${meta.sequencing_run[normal_idx]} --assay ${params.cdm} --id-snp ${params.outdir}/${params.subdir}/QC/${normaljson} "> ${meta.id[normal_idx]}.pairgen
        """

      } else {
        tumor_idx   = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
        tumor_id    = meta.id[tumor_idx]
	    tumorjson    = jsons[tumor_idx]

        """
        echo "--overwrite --sample-id ${meta.id[tumor_idx]} --sequencing-run ${meta.sequencing_run[tumor_idx]} --assay ${params.cdm} --id-snp ${params.outdir}/${params.subdir}/QC/${tumorjson} " > ${meta.id[tumor_idx]}.pairgen
        """
      }
    
    stub:
        def prefix = task.ext.prefix ?: "${meta.id}"
        
        if(meta.id.size() == 2) {
            tumor_idx   = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
            normal_idx  = meta.type.findIndexOf{ it == 'normal' || it == 'N' }
            normal_id   = meta.id[normal_idx]
            tumor_id    = meta.id[tumor_idx]
	        normaljson  = jsons[normal_idx]
	        tumorjson    = jsons[tumor_idx]
             """
             echo "--overwrite --sample-id ${meta.id[tumor_idx]} --sequencing-run ${meta.sequencing_run[tumor_idx]} --assay ${params.cdm} --id-snp ${params.outdir}/${params.subdir}/QC/${tumorjson} " > ${meta.id[tumor_idx]}.pairgen
             
             echo "--overwrite --sample-id ${meta.id[normal_idx]} --sequencing-run ${meta.sequencing_run[normal_idx]} --assay ${params.cdm} --id-snp ${params.outdir}/${params.subdir}/QC/${normaljson} "> ${meta.id[normal_idx]}.pairgen
            """
        } else {
            tumor_idx   = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
            tumor_id    = meta.id[tumor_idx]
	        tumorjson    = jsons[tumor_idx]

            """
            echo "--overwrite --sample-id ${meta.id[tumor_idx]} --sequencing-run ${meta.sequencing_run[tumor_idx]} --assay ${params.cdm} --id-snp ${params.outdir}/${params.subdir}/QC/${tumorjson} " > ${meta.id[tumor_idx]}.pairgen
            """
        }

}

process PROVIDER {
    label 'process_single'
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), file(bam), file(bai)

    
    output:
        tuple val(group), val(meta), file("*.genotypes"),    emit: genotype_checked
        path "versions.yml",                                 emit: versions  // Emit version information in YAML format
    
    when:
        task.ext.when == null || task.ext.when
        
    script:
        def prefix  = task.ext.prefix ?: "${meta.id}"
        def args    = task.ext.args  ?: ""

        """
        provider.pl  --bam $bam  $args  --out $prefix 
    
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
        END_VERSIONS
        """
        
    stub:
        def prefix  = task.ext.prefix ?: "${meta.id}"

        """
        touch ${prefix}.genotypes

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
        END_VERSIONS
        """
}
