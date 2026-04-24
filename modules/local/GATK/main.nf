// SOMATIC CALLING //
process GATKCOV_BAF {
    label 'process_medium_cpus'
    label 'process_high_memory'
    label 'process_more_time'

    input:
        tuple val(group), val(meta), file(cram), file(crai), file(bai)

    output:
        tuple val(group), val(meta), file("*.allelicCounts.tsv"),  emit: gatk_baf
        path "versions.yml",                                       emit: versions
    
    when:
        task.ext.when == null || task.ext.when

    script:
        def args        = task.ext.args   ?: ''
        def prefix      = task.ext.prefix ?: "${meta.id}"
        def avail_mem   = 52288
        if (!task.memory) {
            log.info '[GATK CollectAllelicCounts] Available memory not known - defaulting to 12GB. Specify process memory requirements to change this.'
        } else {
            avail_mem = (task.memory.mega*0.8).intValue()
        }
        """
        export THEANO_FLAGS="base_compiledir=."
        export MKL_NUM_THREADS=${task.cpus}
        export OMP_NUM_THREADS=${task.cpus}
        set +u
        source activate gatk
        gatk --java-options "-Xmx${avail_mem}M" CollectAllelicCounts \\
            -I $cram \\
            $args \\
            -O ${prefix}.allelicCounts.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        export THEANO_FLAGS="base_compiledir=."
        export MKL_NUM_THREADS=${task.cpus}
        export OMP_NUM_THREADS=${task.cpus}
        set +u
        source activate gatk
        touch ${prefix}.allelicCounts.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        END_VERSIONS
        """
}


process GATKCOV_COUNT {
    label 'process_high'
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), file(cram), file(crai), file(bai)

    output:
        tuple val(group), val(meta), file("*.standardizedCR.tsv"), file("*.denoisedCR.tsv"),    emit: gatk_count
        path "versions.yml",                                                                    emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        // Initialize default values
        def args        = task.ext.args      ?: ''
        def args2       = task.ext.args2     ?: ''
        def prefix      = task.ext.prefix    ?: "${meta.id}"
        def sex         = task.ext.sex       ?: "${meta.sex}"    
        def PON         // Initialize the PON variable
        def avail_mem   = 52288

        // Set PON based on sex
        if (sex == "F") {
            PON = params.GATK_PON_FEMALE
        } else {
            PON = params.GATK_PON_MALE
        }
        
        // Handle memory allocation
        if (!task.memory) {
            log.info '[GATK] Available memory not known - defaulting to 50GB. Specify process memory requirements to change this.'
        } else {
            avail_mem = (task.memory.mega * 0.8).intValue()
        } avail_mem = (task.memory.mega*0.8).intValue()
        
        """
        export THEANO_FLAGS="base_compiledir=."
        export MKL_NUM_THREADS=${task.cpus}
        export OMP_NUM_THREADS=${task.cpus}
        set +u
        source activate gatk
        gatk CollectReadCounts \\
            -I ${cram} \\
            $args \\
            -O ${cram}.hdf5

        gatk --java-options "-Xmx${avail_mem}M" DenoiseReadCounts \\
            -I ${cram}.hdf5  \\
            --count-panel-of-normals ${PON} \\
            --standardized-copy-ratios ${prefix}.standardizedCR.tsv \\
            --denoised-copy-ratios ${prefix}.denoisedCR.tsv

        gatk PlotDenoisedCopyRatios \\
            --standardized-copy-ratios ${prefix}.standardizedCR.tsv \\
            --denoised-copy-ratios ${prefix}.denoisedCR.tsv \\
            $args2 \\
            --output . --output-prefix ${prefix}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.id}"
        def sex    = task.ext.sex    ?: "${meta.sex}" 
        """
        export THEANO_FLAGS="base_compiledir=."
        export MKL_NUM_THREADS=${task.cpus}
        export OMP_NUM_THREADS=${task.cpus}
        set +u
        source activate gatk

        echo $sex
        touch ${prefix}.standardizedCR.tsv 
        touch ${prefix}.denoisedCR.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        END_VERSIONS
        """    
}


process GATKCOV_CALL {
    label 'process_high'
    tag "${meta.id[tumor_idx]}"

    input:
        tuple val(group), val(meta), file(allele), file(stdCR), file(denoised)

    output:
        tuple val(group), val(meta), file("*.called.seg"),                  emit: gatkcov_called
        tuple val(group), file("*.modeled.png"),                            emit: gatkcov_plot
        path "versions.yml",                                                emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args   = task.ext.args      ?: ''
        def args2  = task.ext.args2     ?: ''
        def args3  = task.ext.args3     ?: ''

        normalcounts = ""
        tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T'  }
        normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N'  }
        prefix = task.ext.prefix ?: "${meta.id[tumor_idx]}"
        if (meta.id.size() == 2) {
            args = args + " --normal-allelic-counts " + allele[normal_idx]
        }

        def avail_mem   = 52288
        if (!task.memory) {
            log.info '[GATK] Available memory not known - defaulting to 12GB. Specify process memory requirements to change this.'
        } else {
            avail_mem = (task.memory.mega*0.8).intValue()
        }

        """
        export THEANO_FLAGS="base_compiledir=."
        export MKL_NUM_THREADS=${task.cpus}
        export OMP_NUM_THREADS=${task.cpus}
        set +u
        source activate gatk
        gatk --java-options "-Xmx${avail_mem}M" ModelSegments \\
            --denoised-copy-ratios ${denoised[tumor_idx]} \\
            --allelic-counts ${allele[tumor_idx]} \\
            $args \\
            --output . \\
            --output-prefix ${prefix}

        gatk CallCopyRatioSegments \\
            --input ${prefix}.cr.seg \\
            --output ${prefix}.called.seg \\
            $args2

        gatk PlotModeledSegments \\
            --denoised-copy-ratios ${denoised[tumor_idx]} \\
            --allelic-counts ${prefix}.hets.tsv \\
            --segments ${prefix}.modelFinal.seg \\
            $args3 \\
            --output . \\
            --output-prefix ${prefix}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        END_VERSIONS
        """

    stub:
        tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T'  }
        normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N'  }
        prefix = task.ext.prefix ?: "${meta.id[tumor_idx]}"
        """
        export THEANO_FLAGS="base_compiledir=."
        export MKL_NUM_THREADS=${task.cpus}
        export OMP_NUM_THREADS=${task.cpus}
        set +u
        source activate gatk

        touch ${prefix}.called.seg
        touch ${prefix}.modeled.png

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        END_VERSIONS
        """
}


process GATKCOV_CALL_GERMLINE {
    label 'process_high'
    tag "${meta.id[normal_idx]}"

    input:
        tuple val(group), val(meta), file(allele), file(stdCR), file(denoised)

    output:
        tuple val(group), val(meta), file("*.called.seg"),                  emit: gatkcov_germline_called
        tuple val(group), file("*.modeled.png"),                            emit: gatkcov_germlline_plot
        path "versions.yml",                                                emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args   = task.ext.args      ?: ''
        normalcounts = ""
        tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T'  }
        normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N'  }
        prefix = task.ext.prefix ?: "${meta.id[normal_idx]}"
        def avail_mem   = 12288
        if (!task.memory) {
            log.info '[GATK] Available memory not known - defaulting to 12GB. Specify process memory requirements to change this.'
        } else {
            avail_mem = (task.memory.mega*0.8).intValue()
        }

        """
        export THEANO_FLAGS="base_compiledir=."
        export MKL_NUM_THREADS=${task.cpus}
        export OMP_NUM_THREADS=${task.cpus}
        set +u
        source activate gatk
        gatk --java-options "-Xmx${avail_mem}M" ModelSegments \\
            --denoised-copy-ratios ${denoised[normal_idx]} \\
            --output . \\
            --output-prefix ${prefix}

        gatk CallCopyRatioSegments \\
            --input ${prefix}.cr.seg \\
            --output ${prefix}.called.seg

        gatk PlotModeledSegments \\
            --denoised-copy-ratios ${denoised[normal_idx]} \\
            --segments ${prefix}.modelFinal.seg \\
            $args \\
            --output . \\
            --output-prefix ${prefix}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        END_VERSIONS
        """

    stub:
        tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T'  }
        normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N'  }
        prefix = task.ext.prefix ?: "${meta.id[normal_idx]}"
        """
        export THEANO_FLAGS="base_compiledir=."
        export MKL_NUM_THREADS=${task.cpus}
        export OMP_NUM_THREADS=${task.cpus}
        set +u
        source activate gatk

        touch ${prefix}.called.seg
        touch ${prefix}.modeled.png

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        END_VERSIONS
        """
}

