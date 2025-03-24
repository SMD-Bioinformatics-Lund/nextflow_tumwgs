process TNSCOPE {
    label "process_medium"
    tag "$group"

    input:
        tuple val(group), val(meta), file(bams), file(bais), file(bqsr)
        each file(bed)

    output:
        tuple val("tnscope"), val(group), file("tnscope_${bed}.vcf"),   emit: vcfparts_tnscope
        path "versions.yml",                                            emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args     ?: ''
        def args2   = task.ext.args2    ?: ''

        if( meta.id.size() >= 2 ) {
            tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
            normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N' }
            """
            sentieon driver -t ${task.cpus} $args \\
                -i ${bams[tumor_idx]} -q ${bqsr[tumor_idx]} \\
                -i ${bams[normal_idx]} -q ${bqsr[normal_idx]} \\
                --interval $bed --algo TNscope \\
                --tumor_sample ${meta.id[tumor_idx]} --normal_sample ${meta.id[normal_idx]} \\
                $args2 \\
                --min_tumor_allele_frac ${params.tnscope_var_freq_cutoff_p} \\
                tnscope_${bed}.vcf.raw

            filter_tnscope_somatic.pl tnscope_${bed}.vcf.raw ${meta.id[tumor_idx]} ${meta.id[normal_idx]} > tnscope_${bed}.vcf

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
            END_VERSIONS
            """
        }
        else if( meta.id.size() == 1 ) {
            """
            sentieon driver -t ${task.cpus} $args \\
                -i ${bams} -q ${bqsr} \\
                --interval $bed --algo TNscope \\
                --tumor_sample ${meta.id[0]} \\
                $args2 \\
                --min_tumor_allele_frac ${params.tnscope_var_freq_cutoff_up} \\
                tnscope_${bed}.vcf.raw

            filter_tnscope_unpaired.pl tnscope_${bed}.vcf.raw > tnscope_${bed}.vcf

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
            END_VERSIONS
            """ 
        }

    stub:
        """
        touch tnscope_${bed}.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
        END_VERSIONS
        """
}

process BWA_ALIGN_SHARD {
    label 'process_alot'
    tag "${shard}_${meta.id}"

    input:
        tuple val(shard), val(group), val(meta), file(r1), file(r2)

    output:
        tuple val(meta.id), val(group), val(meta), file("${out_bam}"), file("${out_bam}.bai"),  emit: shard_bam
        path "versions.yml",                                                                    emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args     ?: ''                       
        def args2   = task.ext.args2    ?: ''
        def args3   = task.ext.args3    ?: ''

       out_bam = meta.id+"."+shard+"."+meta.type+".bwa.sort.bam"

        if (meta.sub) {
            submbp = params.sample_val / 1000000
            submbp = submbp + "M"
            out_bam = shard+"."+meta.id+"."+meta.type+"."+submbp+".bwa.sort.bam"
        }

        """
        sentieon bwa mem -t ${task.cpus} ${args} '<sentieon fqidx extract ${args2} ${r1} ${r2}' | sentieon util sort $args3 -o ${out_bam} -t ${task.cpus} --sam2bam -i -

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
            bwa: \$(echo \$(sentieon bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
        END_VERSIONS
        """
    
    stub:
        out_bam = shard+"."+meta.id+"."+meta.type+".bwa.sort.bam"

        if (meta.sub) {
            submbp  = params.sample_val / 1000000
            submbp  = submbp + "M"
            out_bam = meta.id+"."+meta.type+"."+submbp+".bwa.sort.bam"
        }

        """
        touch ${out_bam} ${out_bam}.bai

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
            bwa: \$(echo \$(sentieon bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
        END_VERSIONS
        """
}

process BWA_MERGE_SHARDS {
    label 'process_alot'
    label 'scratch'
    label 'stage'
    tag "${meta.id}"
    
    input:
        tuple val(id), val(group), val(meta), file(shard_bams), file(shard_bams_bai)

    output:
        tuple val(group), val(meta), file("${prefix}.bwa.sort.bam"), file("${prefix}.bwa.sort.bam.bai"),  emit: merged_bam
        path "versions.yml",                                        emit: versions
   
    when:
        task.ext.when == null || task.ext.when

    script:
        // group= meta.group.first()
        // meta2 = meta.unique() // should keep the original structure 
        shard_bams = shard_bams.sort(false) { a, b -> a.getBaseName() <=> b.getBaseName() } .join(' ')
        // out_bam = meta.id.first()+"."+meta.type.first()+".bwa.sort.bam"
        prefix = task.ext.prefix ?: "${meta.id}.${meta.type}"

        """
        sentieon util merge -o ${prefix}.bwa.sort.bam ${shard_bams}
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
            bwa: \$(echo \$(sentieon bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
        END_VERSIONS
        """

    stub:
        prefix = task.ext.prefix ?: "${meta.id}.${meta.type}"

        """
        touch ${prefix}.bwa.sort.bam ${prefix}.bwa.sort.bam.bai

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
            bwa: \$(echo \$(sentieon bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
        END_VERSIONS
        """
}

process BAM_CRAM{
    label 'process_alot'
    label 'scratch'
    label 'stage'
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), file(mergedbam), file(mergedbai)

    output:
        tuple val(group), val(meta), file("${prefix}.sort.cram"), file("${prefix}.sort.cram.crai"), file("${prefix}.sort.cram.bai"),    emit: cram_merged
        path "versions.yml",                                                                                                            emit: versions
   

    script:
        def args    = task.ext.args     ?: ''                  
        prefix = task.ext.prefix ?: "${meta.id}.${meta.type}"

        """
        sentieon driver -t ${task.cpus} -i ${mergedbam} ${args} ${prefix}.sort.cram
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
            bwa: \$(echo \$(sentieon bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
        END_VERSIONS
        """
    stub:
        prefix = task.ext.prefix ?: "${meta.id}.${meta.type}"

        """
        touch ${prefix}.sort.cram 
        touch ${prefix}.sort.cram.crai 
        touch ${prefix}.sort.cram.bai

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
            bwa: \$(echo \$(sentieon bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
        END_VERSIONS
        """

}

process MARKDUP {
    label 'process_very_high'
    label 'process_alot'
    label 'scratch'
    label 'stage'
    tag "${meta.id}"
    
    input:
        tuple val(group), val(meta), file(cram), file(crai), file(bai)

    output:
        tuple val(group), val(meta), file("${out_cram}"), file("${out_cram}.crai"),file("${out_cram}.bai"),                             emit: cram_bqsr
        tuple val(group), val(meta), file("*dedup_metrics.txt"),                                                                        emit: cram_metric
        path "versions.yml",                                                                                                            emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args     ?: ""                       
        def args2   = task.ext.args2    ?: ""
        def args3   = task.ext.args3    ?: ""
        def prefix  = task.ext.prefix   ?: ""

        out_cram = meta.id+"."+meta.type+".sort.dedup.cram"

        if (meta.sub) {
            submbp = params.sample_val / 1000000
            submbp = submbp + "M"
            out_cram = meta.id+"."+meta.type+".sort.dedup.cram"

        }
        """
        sentieon driver $args3 -t ${task.cpus} -i $cram --algo LocusCollector $args
        sentieon driver $args3 -t ${task.cpus} -i $cram --algo Dedup $args2 --metrics ${prefix}dedup_metrics.txt $out_cram

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
        END_VERSIONS
        """

    stub:
        def args    = task.ext.args     ?: ""                       
        def args2   = task.ext.args2    ?: ""
        def prefix  = task.ext.prefix   ?: ""

        out_cram = meta.id+"."+meta.type+".dedup.cram"

        if (meta.sub) {
            submbp = params.sample_val / 1000000
            submbp = submbp + "M"
            out_cram = meta.id+"."+meta.type+"."+submbp+".dedup.cram"
        }
        """
        touch ${out_cram}  ${out_cram}.crai ${out_cram}.bai ${prefix}dedup_metrics.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
        END_VERSIONS
        """
}

process REALIGN_INDEL_BQSR {
    label 'process_alot'
    label 'scratch'
    label 'stage'
    tag '${meta.id}'

    input:
        tuple val(group), val(meta), file(cram), file(crai), file(bai)

    output:
        tuple val(group), val(meta), file("${out_cram}"), file("${out_cram}.crai"),file("${out_cram}.bai"),                         emit: cram_bqsr
        tuple val(group), val(meta), file("${out_cram}"), file("${out_cram}.crai"),file("${out_cram}.bai"), file("*.bqsr.table"),   emit: cram_varcall
        path "versions.yml",                                                                    emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args     ?: ""   // reference and common arguments for driver
        def args2   = task.ext.args2    ?: ""   // algo specific arguments
        def args3   = task.ext.args3    ?: ""   // algo specific arguments
        def prefix  = task.ext.prefix   ?: "${meta.id}"
        out_cram = meta.id+"."+meta.type+".sort.dedup.realign.cram"

        """
        sentieon driver $args -t ${task.cpus} -i $cram --algo Realigner $args2 $out_cram
        sentieon driver $args -t ${task.cpus} -i $out_cram --algo QualCal $args3 ${prefix}.bqsr.table

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
        END_VERSIONS
        """

    stub:
        out_cram = meta.id+"."+meta.type+".sort.dedup.realign.cram"
        """
        touch ${out_cram} ${out_cram}.crai ${out_cram}.bai 
        touch ${out_cram}.bqsr.table

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
        END_VERSIONS
        """
}

process SENTIEON_QC {
    label 'process_medium'
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), file(cram), file(crai), file(bai), file(dedup)

    output:
        tuple val(group), val(meta), file("*_${meta.type}.QC"),                        emit: qc_cdm
        path "versions.yml",                                                           emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args     ?: ""   // reference 
        def prefix  = task.ext.prefix   ?: "${meta.id}"
        """
        sentieon driver $args \\
            -t ${task.cpus} -i ${cram} \\
            --algo MeanQualityByCycle mq_metrics.txt --algo QualDistribution qd_metrics.txt \\
            --algo GCBias --summary gc_summary.txt gc_metrics.txt --algo AlignmentStat aln_metrics.txt \\
            --algo InsertSizeMetricAlgo is_metrics.txt \\
            --algo WgsMetricsAlgo wgs_metrics
            
        qc_sentieon.pl ${meta.id}_${meta.type} wgs > ${prefix}_${meta.type}.QC
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
        END_VERSIONS
        """

    stub:
        def prefix  = task.ext.prefix   ?: "${meta.id}"
        """
        touch ${prefix}_${meta.type}.QC

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
        END_VERSIONS
        """
}

process CRAM_TO_BAM {
    label 'process_alot'
    label 'scratch'
    label 'stage'
    tag '${meta.id}'

    input:
        tuple val(group), val(meta), file(cram), file(crai), file(bai)

    output:
        tuple val(group), val(meta), file("${out_bam}"), file("${out_bam}.bai"), emit: bam_bqsr
        path "versions.yml",                                                     emit: versions

    script:
        def args    = task.ext.args ?: "" 
        out_bam = "${meta.id}.${meta.type}.sort.dedup.realign.bam"

        """
        samtools view $args -@ ${task.cpus} -o ${out_bam}  ${cram}
        samtools index -@ ${task.cpus} ${out_bam}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(echo \$(samtools 2>&1) | sed 's/.*Version: //; s/ .*//')
        END_VERSIONS
        """
    
    stub:
        out_bam = meta.id+"."+meta.type+".sort.dedup.realign.bam"

        """
        touch ${out_bam} ${out_bam}.bai

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(echo \$(samtools 2>&1) | sed 's/.*Version: //; s/ .*//')
        END_VERSIONS
        """

}

process DNASCOPE {
    label 'process_alot'
    label 'scratch'
    label 'stage'
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), file(cram), file(crai), file(bai)

    output:
        tuple val(group), val(meta), file("${out_vcf}"), file("${out_vcf}.tbi"),        emit: dnascope_vcf
        path "versions.yml",                                                            emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args     ?: ""   // reference and common arguments for driver
        def args2   = task.ext.args2    ?: ""   // algo specific arguments
        def args3   = task.ext.args3    ?: ""   // algo specific arguments
        def prefix  = task.ext.prefix   ?: "${meta.id}"
        out_vcf     = meta.id+"."+meta.type+".dnascope.vcf.gz"

        """
        sentieon driver -t ${task.cpus} $args -i $cram $args2 --algo DNAscope $args3 $out_vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
        END_VERSIONS
        """

    stub:
        out_vcf = meta.id+"."+meta.type+".dnascope.vcf.gz"
        """
        touch ${out_vcf} ${out_vcf}.tbi

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
        END_VERSIONS
        """
}