process CONCATENATE_VCFS {
    label "process_single"
    tag "$group"

    input:
        tuple val(vc), val(group), file(vcfs)

    output:
        tuple val(group), val(vc), file("*_${vc}.vcf.gz"),    emit: concatenated_vcfs
        path "versions.yml",                                  emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: "${group}"
        """
        vcf-concat $vcfs | vcf-sort -c | gzip -c > ${vc}.concat.vcf.gz
        vt decompose ${vc}.concat.vcf.gz -o ${vc}.decomposed.vcf.gz
        vt index ${vc}.decomposed.vcf.gz
        vt sort -m chrom ${vc}.decomposed.vcf.gz -o ${vc}.decomposed.sorted.vcf.gz
        vt normalize ${vc}.decomposed.sorted.vcf.gz $args | vt uniq - -o ${prefix}_${vc}.vcf.gz
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            vcftools: \$(echo \$(vcftools --version 2>&1) | sed 's/^.*VCFtools (//;s/).*//')
            vt-decompose: \$(echo \$(vt decompose 2>&1) | sed 's/.*decompose v//; s/ .*//')
            vt-normalize: \$(echo \$(vt normalize 2>&1) | sed 's/.*normalize v//; s/ .*//')
        END_VERSIONS
        """
        
    stub:
        def prefix = task.ext.prefix ?: "${group}"
        """
        touch ${prefix}_${vc}.vcf.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            vcftools: \$(echo \$(vcftools --version 2>&1) | sed 's/^.*VCFtools (//;s/).*//')
            vt-decompose: \$(echo \$(vt decompose 2>&1) | sed 's/.*decompose v//; s/ .*//')
            vt-normalize: \$(echo \$(vt normalize 2>&1) | sed 's/.*normalize v//; s/ .*//')
        END_VERSIONS
        """
}


process AGGREGATE_VCFS {
    label "process_single"
    tag "$group"

    input:
        tuple val(group), val(vc), file(vcfs), val(meta)

    output:
        tuple val(group), val(meta), file("*.agg.vcf"), emit: vcf_concat
        path "versions.yml",                            emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix = task.ext.prefix ?: "${group}"

        sample_order = meta.id[0]
        if( meta.id.size() >= 2 ) {
            tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
            normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N' }
            sample_order = meta.id[tumor_idx]+","+meta.id[normal_idx]
        }

        """
        aggregate_vcf.pl --vcf ${vcfs.sort(false) { a, b -> a.getBaseName() <=> b.getBaseName() }.join(",")} --sample-order ${sample_order} |vcf-sort -c > ${prefix}.agg.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${group}"

        if( meta.id.size() >= 2 ) {
            tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
            normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N' }
            sample_order = meta.id[tumor_idx]+","+meta.id[normal_idx]
            """
            echo tumor:${meta.id[tumor_idx]} normal:${meta.id[normal_idx]}
            touch ${prefix}.agg.vcf

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
            END_VERSIONS
            """
        }
        else {
            """
            touch ${group}.agg.vcf

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
            END_VERSIONS
            """
        }
}
