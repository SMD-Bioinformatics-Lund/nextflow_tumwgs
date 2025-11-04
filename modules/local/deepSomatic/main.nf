process DEEPSOMATIC {
    label 'process_alot'
    label 'scratch'
    label 'stage'
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), file(bams), file(bais)

    output:
        tuple val(group), val("deepSomatic"), file("*.vcf.gz"), emit: deepSomatic_vcf
        path "versions.yml", emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: "${group}"

        if (meta.id.size() >= 2) {
            def tumor_idx  = meta.type.findIndexOf { it == 'tumor' || it == 'T' }
            def normal_idx = meta.type.findIndexOf { it == 'normal' || it == 'N' }

            def output_vcf = "${meta.id[tumor_idx]}_deepsomatic_output.vcf.gz"
            def log_dir    = "${meta.id[tumor_idx]}/logs"
            def interm_dir = "${meta.id[tumor_idx]}/intermediate_results_dir"

            """
            run_deepsomatic \\
                --model_type=WGS \\
                $args \\
                --reads_normal=${bams[normal_idx]} \\
                --reads_tumor=${bams[tumor_idx]} \\
                --output_vcf=${output_vcf} \\
                --sample_name_tumor="${meta.id[tumor_idx]}" \\
                --sample_name_normal="${meta.id[normal_idx]}" \\
                --num_shards=${task.cpus} \\
                --logging_dir=${log_dir} \\
                --intermediate_results_dir=${interm_dir}

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                DeepSomatic: \$(run_deepsomatic --version 2>/dev/null |sed -e "s/DeepSomatic: //g")
            END_VERSIONS
            """
        }
        else if (meta.id.size() == 1) {
            def output_vcf = "${meta.id[0]}_deepsomatic_output.vcf.gz"
            def log_dir    = "${meta.id[0]}/logs"
            def interm_dir = "${meta.id[0]}/intermediate_results_dir"

            """
            run_deepsomatic \\
                --model_type=WGS_TUMOR_ONLY \\
                $args \\
                --reads_tumor=${bams[0]} \\
                --output_vcf=${output_vcf} \\
                --sample_name_tumor="${meta.id[0]}" \\
                --num_shards=${task.cpus} \\
                --logging_dir=${log_dir} \\
                --intermediate_results_dir=${interm_dir} \\
                --use_default_pon_filtering=true

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                DeepSomatic: \$(run_deepsomatic --version 2>&1 | sed -e "s/: DeepVariant version //g")
            END_VERSIONS
            """
        }

    stub:
        if (meta.id.size() >= 2) {
            def tumor_idx  = meta.type.findIndexOf { it == 'tumor' || it == 'T' }
            def normal_idx = meta.type.findIndexOf { it == 'normal' || it == 'N' }

			"""

            echo "Group: $group" > stub.log
            echo "Tumor sample: ${meta.id[tumor_idx]}" >> stub.log
            echo "Normal sample: ${meta.id[normal_idx]}" >> stub.log
            echo "Tumor BAM: ${bams[tumor_idx]}" >>  stub.log
            echo "Normal BAM: ${bams[normal_idx]}" >> stub.log
            echo "Arguments: $args" >> stub.log

            echo "${meta.id[tumor_idx]}_deepsomatic_output.vcf.gz" > ${meta.id[tumor_idx]}_deepsomatic_output.vcf.gz

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                DeepSomatic: \$(run_deepsomatic --version 2>/dev/null |sed -e "s/DeepSomatic: //g")
            END_VERSIONS
			"""
        }
        else {
			"""
            echo "Group: $group" > stub.log
            echo "Tumor sample: ${meta.id[0]}" >> stub.log
            echo "Tumor BAM: ${bams[0]}" >> stub.log
            echo "Arguments: $args" >> stub.log

            echo "${meta.id[0]}_deepsomatic_output.vcf.gz" > ${meta.id[0]}_deepsomatic_output.vcf.gz

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                DeepSomatic: stub
            END_VERSIONS
			"""
        }
}
