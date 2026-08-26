process GERMLINE_EVALUATE {
    label "process_single"
    tag "$group"

    input:
        tuple val(group), val(meta), file(vcf)

    output:
        tuple val(group), val(meta), file("*.germline.ranked.vcf"),   emit: vcf_ranked
        path "versions.yml",                                          emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args   = task.ext.args   ?: ''
        def prefix = task.ext.prefix ?: "${vcf.baseName}"
        tumor_idx  = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }

        """
        germline_evaluate.pl --vcf $vcf --tumor-id ${meta.id[tumor_idx]} $args > ${prefix}.germline.ranked.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${group}"
        tumor_idx  = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
        """
        echo ${meta.id[tumor_idx]}
        touch ${prefix}.germline.ranked.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
        END_VERSIONS
        """
}
