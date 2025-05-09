process FASTP {
    label 'process_very_high'
    label 'max_errors'
    label 'scratch'
    label 'stage'
    tag "${meta.id}"


    input:
        tuple val(group), val(meta), file(r1), file(r2)

    output:
        tuple val(group), val(meta), file("*.cleaned_R1.fastq.gz"), file("*.cleaned_R2.fastq.gz"),  emit: fastq_trimmed
        tuple val(group), file("*.cleaned.fastp.json"), file("*.cleaned.html"),                     emit: fastq_stats
        path "versions.yml",                                                                        emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: "${meta.id}.${meta.type}"
        """
        fastp -w ${task.cpus} \\
        $args \\
        -i $r1 -I $r2 \\
        -o ${prefix}.cleaned_R1.fastq.gz \\
        -O ${prefix}.cleaned_R2.fastq.gz \\
        -h ${prefix}.cleaned.html \\
        -j ${prefix}.cleaned.fastp.json

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.id}.${meta.type}"
        """
        touch ${prefix}.cleaned_R1.fastq.gz
        touch ${prefix}.cleaned_R2.fastq.gz
        touch ${prefix}.cleaned.html
        touch ${prefix}.cleaned.fastp.json

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
        END_VERSIONS
        """
}