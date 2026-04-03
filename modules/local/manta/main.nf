process MANTA {
    label 'process_alot'
    tag "$group"
    
    input:
        tuple val(group), val(meta), file(cram), file(crai), file(bai)

    output:
        tuple val(group),val(meta), file("${prefix}_manta.vcf"),                emit: manta_vcf_tumor
        tuple val(group),val(meta), file("${prefix2}_manta.vcf"),               emit: manta_vcf_normal
        path "versions.yml",                                                    emit: versions
        

    when:
        task.ext.when == null || task.ext.when
    
    script:
        def args    = task.ext.args  ?: ""
        def args2   = task.ext.args2 ?: ""
        tumor_idx   = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
        normal_idx  = meta.type.findIndexOf{ it == 'normal' || it == 'N' }
        normal      = cram[normal_idx]
        normal_id   = meta.id[normal_idx]
        tumor       = cram[tumor_idx]
        tumor_id    = meta.id[tumor_idx]
        prefix      = task.ext.prefix  ?: tumor_id
        prefix2     = task.ext.prefix2 ?: normal_id

        if(meta.id.size() == 2) { 
            """
            configManta.py \\
                --tumorBam $tumor \\
                --normalBam $normal \\
                $args \\
                --runDir .
            python runWorkflow.py $args2
            mv ./results/variants/somaticSV.vcf.gz ${prefix}_manta.${tumor}.vcf.gz
            mv ./results/variants/diploidSV.vcf.gz ${prefix2}_manta.${normal}.vcf.gz
            gunzip ${prefix}_manta.${tumor}.vcf.gz
            gunzip ${prefix2}_manta.${normal}.vcf.gz
            mv ${prefix}_manta.${tumor}.vcf ${prefix}_manta.vcf
            mv ${prefix2}_manta.${normal}.vcf ${prefix2}_manta.vcf


            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                manta: \$(configManta.py --version)
                python: \$(python --version 2>&1 | sed -e 's/Python //g')
            END_VERSIONS
            """
        }
        else {
            """
            set +eu
            source activate py2
            set -eu
            configManta.py \\
                --tumorBam $cram \\
                $args \\
                --runDir .
            python runWorkflow.py $args2
            mv results/variants/tumorSV.vcf.gz ${prefix}_manta.vcf.gz
            gunzip ${prefix}_manta.vcf.gz

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                manta: \$(configManta.py --version)
                python: \$(python --version 2>&1 | sed -e 's/Python //g')
            END_VERSIONS
            """
        }

    stub:
        tumor_idx   = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
        normal_idx  = meta.type.findIndexOf{ it == 'normal' || it == 'N' }
        normal      = cram[normal_idx]
        normal_id   = meta.id[normal_idx]
        tumor       = cram[tumor_idx]
        tumor_id    = meta.id[tumor_idx]
        prefix      = task.ext.prefix  ?: tumor_id
        prefix2     = task.ext.prefix2 ?: normal_id
        if(meta.id.size() == 2) {
            """
            set +eu
            source activate py2
            set -eu
            touch ${prefix}_manta.vcf 
            touch ${prefix2}_manta.vcf

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                manta: \$(configManta.py --version)
                python: \$(python --version 2>&1 | sed -e 's/Python //g')
            END_VERSIONS
            """
        }
        else {
            """
            set +eu
            source activate py2
            set -eu
            touch ${prefix}_manta.vcf 

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                manta: \$(configManta.py --version)
                python: \$(python --version 2>&1 | sed -e 's/Python //g')
            END_VERSIONS
            """
        }

}


process MANTA_SV {
    label 'process_medium'
    tag "${meta.id}"
    
    input:
        tuple val(group), val(meta), file(vcf)

    output:
        tuple val(group),val(meta), file("*.manta.SV.pass.vcf"),                emit: inv_vcf_tumor
        path "versions.yml",                                                    emit: versions
        
    when:
        task.ext.when == null || task.ext.when

    script:
        def args        = task.ext.args ?: ''
        def prefix      = task.ext.prefix ?: "${group}"
    
        """
        ## convertInversion.py $args ${vcf} > ${prefix}.manta.all.vcf # need to rethink the inversion filtering, as the current one is missed by snpeff annotation
        grep -e \$'\\tPASS\\t' -e '^#' ${vcf} |grep -Ev 'GL000|hs37d5' > ${prefix}.manta.SV.pass.vcf
    
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            manta: \$(configManta.py --version)
            python: \$(python --version 2>&1 | sed -e 's/Python //g')
        END_VERSIONS    
        """
    stub:
        def args        = task.ext.args ?: ''
        def prefix      = task.ext.prefix ?: "${group}" 
        """
        touch ${prefix}.manta.SV.pass.vcf
        echo ${args}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            manta: \$(configManta.py --version)
            python: \$(python --version 2>&1 | sed -e 's/Python //g')
        END_VERSIONS 
        """
}




