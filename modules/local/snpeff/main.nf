process SNPEFF {
    label 'process_medium'
    tag "$group"

    input:
        tuple val(group), val(meta), file(vcf)

    output:
        tuple val(group), val(meta), file("*.merged.annotated.vcf"),            emit: snpeff_vcf
        path "versions.yml",                                                    emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args        = task.ext.args ?: ''
        def prefix      = task.ext.prefix ?: "${group}"
        def avail_mem   = 6144
        if (!task.memory) {
            log.info '[snpEff] Available memory not known - defaulting to 6GB. Specify process memory requirements to change this.'
        } else {
            avail_mem = (task.memory.mega*0.8).intValue()
        }
        """
        snpEff -Xmx${avail_mem}M $args ${vcf} > ${prefix}.merged.annotated.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            snpEff: \$(echo \$(snpEff -version 2>&1) | grep 'SnpEff ' | sed 's/.*SnpEff //; s/ .*\$//')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${group}"
        """
        touch ${prefix}.merged.annotated.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            snpEff: \$(echo \$(snpEff -version 2>&1) | grep 'SnpEff ' | sed 's/.*SnpEff //; s/ .*\$//')
        END_VERSIONS
        """ 
}


process SNPEFF_SV_ANN {
    label 'process_medium'
    tag "$group"

    input:
        tuple val(group), val(meta), file(vcf)

    output:
        tuple val(group), val(meta), file("*.BND.annotated.vcf"), file("*.TANDEM.INS.DEL.SV_annotated.vcf"),    emit: snpeff_BND_TANDEM
        path "versions.yml",                                                                                    emit: versions
    
    when:
        task.ext.when == null || task.ext.when

    script:
        def args        = task.ext.args ?: ''
        def prefix      = task.ext.prefix ?: "${group}"
        def avail_mem   = 6144
        if (!task.memory) {
            log.info '[snpEff] Available memory not known - defaulting to 6GB. Specify process memory requirements to change this.'
        } else {
            avail_mem = (task.memory.mega*0.8).intValue()
        }
        """
        snpEff -Xmx${avail_mem}M $args ${vcf} > ${prefix}.SV.annotated.vcf
        grep -e '^#' -e 'MantaBND:' ${prefix}.SV.annotated.vcf >  ${prefix}.BND.annotated.vcf
        grep -v 'MantaBND:' ${prefix}.SV.annotated.vcf | grep -e '^#' -e 'MantaINV:' >  ${prefix}.INV.annotated.vcf
        grep -v 'MantaBND:' ${prefix}.SV.annotated.vcf | grep -v 'MantaINV:'  >  ${prefix}.TANDEM.INS.DEL.SV_annotated.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            snpEff: \$(echo \$(snpEff -version 2>&1) | grep 'SnpEff ' | sed 's/.*SnpEff //; s/ .*\$//')
        END_VERSIONS
        """

    stub:
        def args        = task.ext.args ?: ''
        def prefix      = task.ext.prefix ?: "${group}"
        """
        touch ${prefix}.BND.annotated.vcf
        touch ${prefix}.TANDEM.INS.DEL.SV_annotated.vcf


        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            snpEff: \$(echo \$(snpEff -version 2>&1) | grep 'SnpEff ' | sed 's/.*SnpEff //; s/ .*\$//')
        END_VERSIONS
        """ 
}


process SNPEFF_SV_FINAL {
    label 'process_medium'
    tag "$group"

    input:
        tuple val(group), val(meta), file(vcf)

    output:
        tuple val(group), val(meta), file("*.BND.annotated.vcf"),file("*.TANDEM.SV_annotated.vcf"), emit: snpeff_BND_TANDEM
        path ("CMD-fusion.annotated.txt"),                             emit: snpeff_CMD
        path "versions.yml",                                                                        emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args        = task.ext.args ?: ''
        def prefix      = task.ext.prefix ?: "${group}"
        def avail_mem   = 6144
        if (!task.memory) {
            log.info '[snpEff] Available memory not known - defaulting to 6GB. Specify process memory requirements to change this.'
        } else {
            avail_mem = (task.memory.mega*0.8).intValue()
        }
        """
        snpEff -Xmx${avail_mem}M $args ${vcf} > ${prefix}.SV.annotated.vcf
        grep -e '^#' -e 'MantaBND:' ${prefix}.SV.annotated.vcf >  ${prefix}.BND.annotated.vcf
        grep -v 'MantaBND:' ${prefix}.SV.annotated.vcf  | grep -e '^#' -e 'MantaINV:' >  ${prefix}.INV.annotated.vcf
        grep -v 'MantaBND:' ${prefix}.SV.annotated.vcf | grep -v 'MantaINV:' | grep -e '^#' -e 'TANDEM' >  ${prefix}.TANDEM.SV_annotated.vcf

        grep -ve '^#' ${prefix}.BND.annotated.vcf |grep CMD > CMD-fusion.txt
        grep -ve '^#' ${prefix}.TANDEM.SV_annotated.vcf |grep BRAF |less -Sx20 >> CMD-fusion.txt
        
        awk -F '\t' '{
            OFS="\\t";
            print \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$(NF-3)";PANEL=fusion|somatic|both", \$9, \$10, \$11
        }' CMD-fusion.txt > CMD-fusion.annotated.txt


        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            snpEff: \$(echo \$(snpEff -version 2>&1) | grep 'SnpEff ' | sed 's/.*SnpEff //; s/ .*\$//')
        END_VERSIONS
        """

    stub:
        def args        = task.ext.args ?: ''
        def prefix      = task.ext.prefix ?: "${group}"
        """
        echo ${args}
        touch ${prefix}.BND.annotated.vcf
        touch ${prefix}.TANDEM.SV_annotated.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            snpEff: \$(echo \$(snpEff -version 2>&1) | grep 'SnpEff ' | sed 's/.*SnpEff //; s/ .*\$//')
        END_VERSIONS
        """         
}


process COMBINE_FUSIONS { 
    label 'process_medium'
    tag "$group"

    input:
        tuple val(group), val(meta), file(vcf)
        tuple val(group), val(meta), file(bnd),file(tandem)

    output:
        tuple val(group), file("*.final.fusions.vcf"),                          emit: sv_CMD
        path "versions.yml",                                                    emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${group}"
        """
        grep -ve '^#' ${bnd} | grep 'CMD' > a.txt || true
        grep -ve '^#' ${tandem} | grep -f ${params.FUSIONS_CNV} > b.txt || true
        cat a.txt b.txt > CMD_fusion.txt

        if [[ -s CMD_fusion.txt ]]; then
                awk -F "\\t" 'BEGIN { OFS="\\t" } { print \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$(NF-3) ";PANEL=fusion|somatic|both", \$9, \$10, \$11 }' CMD_fusion.txt > Selected.txt
        else
            touch Selected.txt
        fi

        cat ${vcf} Selected.txt > a.vcf        
        echo '##INFO=<ID=PANEL,Number=1,Type=String,Description="Panel origin of the variant">' > extra_header.txt
        bcftools annotate -h extra_header.txt a.vcf ${args} fixed.vcf
        bcftools sort fixed.vcf -o ${prefix}.final.fusions.vcf	

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bcftools: \$(echo \$(bcftools --version 2>&1) | sed 's/bcftools //; s/ .*//')
        END_VERSIONS 
        """

    stub:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${group}"
        """
        echo '##INFO=<ID=PANEL,Number=1,Type=String,Description="Panel origin of the variant">' > extra_header.txt
        cat extra_header.txt
        touch ${prefix}.final.fusions.vcf	

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bcftools: \$(echo \$(bcftools --version 2>&1) | sed 's/bcftools //; s/ .*//')
        END_VERSIONS
        """
}