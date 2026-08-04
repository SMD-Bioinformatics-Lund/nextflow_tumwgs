process PINDEL_CONFIG {
    label "process_single"
    tag "$group"

    input:
        tuple val(group), val(meta), file(cram), file(crai), file(bai), file(ins_size)

    output:
        tuple val(group), val(meta), file("*.pindel_config"),   emit: pindel_config

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix	= task.ext.prefix	?:  "${group}" 

        if( meta.id.size() >= 2 ) {
            tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
            normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N' }
            ins_tumor = ins_size[tumor_idx]
            ins_normal = ins_size[normal_idx]
            cram_tumor = cram[tumor_idx]
            cram_normal = cram[normal_idx]
            id_tumor = meta[tumor_idx].id
            id_normal = meta[normal_idx].id

            """
            INS_T="\$(sed -n '3p' $ins_tumor | cut -f 1 | awk '{print int(\$1+0.5)}')"
            INS_N="\$(sed -n '3p' $ins_normal | cut -f 1 | awk '{print int(\$1+0.5)}')"
            echo "$cram_tumor\t\$INS_T\t$id_tumor" > ${prefix}.pindel_config
            echo "$cram_normal\t\$INS_N\t$id_normal" >> ${prefix}.pindel_config
            """
        }
        else {
            tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
            ins_tumor = ins_size[tumor_idx]
            cram_tumor = cram[tumor_idx]
            id_tumor = meta[tumor_idx].id

            """
            INS_T="\$(sed -n '3p' $ins_tumor | cut -f 1 | awk '{print int(\$1+0.5)}')"
            echo "$cram_tumor\t\$INS_T\t$id_tumor" > ${prefix}.pindel_config
            """
        }

    stub:
        def prefix	= task.ext.prefix	?:  "${group}"
        """
        touch ${prefix}.pindel_config
        """
}

process PINDEL_CALLING {
    label "process_medium"
    tag "$group"

    input:
        tuple val(group), val(meta), file(cram), file(crai), file(bai), file(ins_size)
        tuple val(group), val(meta), file(pindel_config)

    output:
        tuple val(group), val("pindel"), file("*_pindel.vcf"),  emit: pindel_vcf
        path "versions.yml",                                    emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args     ?: ''
        def args2   = task.ext.args2	?: ''
        def prefix	= task.ext.prefix	?:  "${group}" 
        """
        pindel $args -i $pindel_config -o tmpout -T ${task.cpus}
        pindel2vcf $args2 -P tmpout -v ${prefix}_pindel_unfilt.vcf
        filter_pindel_somatic.pl ${prefix}_pindel_unfilt.vcf ${prefix}_pindel.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            pindel: \$(echo \$(pindel --version 2>&1)  | grep 'Pindel version' | sed 's/.*Pindel version //' | sed 's/, .*//g')
        END_VERSIONS
        """

    stub:
        def prefix	= task.ext.prefix	?:  "${group}" 
        """
        touch ${prefix}_pindel.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            pindel: \$(echo \$(pindel --version 2>&1)  | grep 'Pindel version' | sed 's/.*Pindel version //' | sed 's/, .*//g')
        END_VERSIONS
        """
}