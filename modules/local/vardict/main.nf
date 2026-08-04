process VARDICT {
    label "process_low"
    tag "$group"

    input:
        tuple val(group), val(meta), file(bams), file(bais)
        each file(bed)

    output:
        tuple val("vardict"), val(group), file("vardict_${bed}.vcf"),   emit: vcfparts_vardict
        path "versions.yml",                                            emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''

        if( meta.id.size() >= 2 ) {

            tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
            normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N' }

            """
            vardict-java -N ${meta.id[tumor_idx]} \\
            -b "${bams[tumor_idx]}|${bams[normal_idx]}" \\
            -U $bed \\
            -f $params.vardict_var_freq_cutoff_p \\
            $args \\
            | testsomatic.R \\
            | var2vcf_paired.pl -N "${meta.id[tumor_idx]}|${meta.id[normal_idx]}" \\
            -f $params.vardict_var_freq_cutoff_p > vardict_${bed}.vcf.raw

            filter_vardict_somatic_wgs.pl vardict_${bed}.vcf.raw ${meta.id[tumor_idx]} ${meta.id[normal_idx]} | grep -v 'FAIL_' > vardict_${bed}.vcf

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                vardict: \$( realpath \$( command -v vardict-java ) | sed 's/.*java-//;s/-.*//' )
                perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
            END_VERSIONS
            """
        }
        else if( meta.id.size() == 1 ) {
            """
            vardict-java -N ${meta.id[0]} \\
            -b ${bams[0]} \\
            -U $bed \\
            -f $params.vardict_var_freq_cutoff_up \\
            $args \\
            | teststrandbias.R \\
            | var2vcf_valid.pl -N ${meta.id[0]} \\
            -E -f 0.01 > vardict_${bed}.vcf.raw

            filter_vardict_unpaired.pl vardict_${bed}.vcf.raw > vardict_${bed}.vcf

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                vardict: \$( realpath \$( command -v vardict-java ) | sed 's/.*java-//;s/-.*//' )
                perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
            END_VERSIONS
            """
        }

    stub:
        if( meta.id.size() >= 2 ) {
            tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
            normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N' }
            """
            echo tumor:${bams[tumor_idx]} ${meta.id[tumor_idx]} normal:${bams[normal_idx]} ${meta.id[normal_idx]}
            touch vardict_${bed}.vcf

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                vardict: \$( realpath \$( command -v vardict-java ) | sed 's/.*java-//;s/-.*//' )
                perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
            END_VERSIONS
            """
        }
        else {
            """
            echo tumor:$bams
            touch vardict_${bed}.vcf

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                vardict: \$( realpath \$( command -v vardict-java ) | sed 's/.*java-//;s/-.*//' )
                perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
            END_VERSIONS
            """
        }
}