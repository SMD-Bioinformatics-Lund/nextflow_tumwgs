process SEQTK {
    label 'scratch'
    label 'stage'
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), file(r1), file(r2)

    output:
        tuple val(group), val(meta), file("*_R1_subsample_${meta.sub}.fastq.gz"), file("*_R2_subsample_${meta.sub}.fastq.gz"),  emit: fastq_sub
        path "versions.yml",                                                                                                    emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args ?: ''
        def prefix  = task.ext.prefix ?: "${meta.clarity_sample_id}_${meta.clarity_pool_id}"
        """
        seqtk sample $args $r1 ${meta.sub} | gzip --no-name > ${prefix}_R1_subsample_${meta.sub}.fastq.gz &
        seqtk sample $args $r2 ${meta.sub} | gzip --no-name > ${prefix}_R2_subsample_${meta.sub}.fastq.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            seqtk: \$(echo \$(seqtk 2>&1) | sed 's/.*Version: //; s/ .*//')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.clarity_sample_id}_${meta.clarity_pool_id}"
        """
        touch ${prefix}_R1_subsample_${meta.sub}.fastq.gz
        touch ${prefix}_R2_subsample_${meta.sub}.fastq.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            seqtk: \$(echo \$(seqtk 2>&1) | sed 's/.*Version: //; s/ .*//')
        END_VERSIONS
        """
}