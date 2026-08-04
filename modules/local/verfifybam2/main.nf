process VERIFYBAMID2 {
    label 'process_medium'
    tag "${meta.id}"
    label 'stage'

	input:
		tuple val(group), val(meta), file(cram), file(crai), file(bai)

    output:
        tuple val(group), val(meta), file("*.contamination.json"),            emit: contamination_json
        tuple val(group), file("*.result.Ancestry"), file("*.result.selfSM"), emit: results
        path "versions.yml",                                                  emit: versions

    script:
        def args    = task.ext.args     ?: ""   // reference 
        def args2   = task.ext.args2    ?: ""   // loci to check
        def prefix  = task.ext.prefix   ?: "${meta.id}"
        """
        verifybamid2 \
            $args \
            $args2 \
            --BamFile $cram
        
        mv result.selfSM ${prefix}.result.selfSM
        mv result.Ancestry ${prefix}.result.Ancestry
        tail -n +2 ${prefix}.result.selfSM |cut -f 7 > ${prefix}.contamination.value
        value=\$(cat ${prefix}.contamination.value)
        echo "{ \\"contamination\\": \\"\$value\\" }" > ${prefix}.contamination.json

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            verifybamid: \$(echo \$(verifybamid2 -h 2>&1 | grep Version | sed "s/Version://"))
        END_VERSIONS
        """

    stub:
        def prefix  = task.ext.prefix   ?: "${meta.id}"
        """
        touch ${prefix}.contamination.json
        touch  ${prefix}.result.Ancestry ${prefix}.result.selfSM

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            verifybamid: \$(echo \$(verifybamid2 -h 2>&1 | grep Version | sed "s/Version://"))
        END_VERSIONS
        """
}