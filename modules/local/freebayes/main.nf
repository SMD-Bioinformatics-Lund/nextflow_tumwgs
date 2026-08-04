process FREEBAYES {
    label "process_low"
    tag "$group"
    
    input:
        tuple val(group), val(meta), file(bams), file(bais)
        each file(bed)

    output:
        tuple val("freebayes"),  val(group), file("freebayes_${bed}.vcf"),  emit: vcfparts_freebayes
        path "versions.yml",                                                emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args                ?: ''
        def args2   = task.ext.args2               ?: ''
        def args3   = task.ext.args3               ?: ''

        if( meta.id.size() >= 2 ) {

            tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
            normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N' }

            """
            freebayes -t $bed \\
            $args \\
            -F ${params.fb_var_freq_cutoff_p} \\
            ${bams[tumor_idx]} \\
            ${bams[normal_idx]} > freebayes_${bed}.vcf.raw

            # vcffilter $args2 freebayes_${bed}.vcf.raw \\
            # | vcffilter $args3 \\
            # | vcfglxgt > freebayes_${bed}.filt1.vcf

            filter_freebayes_somatic_wgs.pl freebayes_${bed}.vcf.raw ${meta.id[tumor_idx]} ${meta.id[normal_idx]} |grep -v 'FAIL_'  > freebayes_${bed}.vcf

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                freebayes: \$(echo \$(freebayes --version 2>&1) | sed 's/version:\s*v//g' )
                vcffilter: \$(echo \$( vcffilter -h 2>&1) | grep 'vcflib' | sed 's/ filter.*\$//g' | sed  's/.* //g' )
                perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
            END_VERSIONS
            """
        }
        else if( meta.id.size() == 1 ) {
            """
            freebayes -t $bed \\
            $args \\
            -F ${params.fb_var_freq_cutoff_up} \\
            $bams > freebayes_${bed}.vcf.raw

            # vcffilter $args2 freebayes_${bed}.vcf.raw \\
            # | vcffilter $args3 \\
            # | vcfglxgt > freebayes_${bed}.filt1.vcf

            filter_freebayes_unpaired.pl freebayes_${bed}.vcf.raw |grep -v 'FAIL_' > freebayes_${bed}.vcf

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                freebayes: \$(echo \$(freebayes --version 2>&1) | sed 's/version:\s*v//g' )
                vcffilter: \$(echo \$( vcffilter -h 2>&1) | grep 'vcflib' | sed 's/ filter.*\$//g' | sed 's/.* //g' )
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
            touch freebayes_${bed}.vcf

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                freebayes: \$(echo \$(freebayes --version 2>&1) | sed 's/version:\s*v//g' )
                vcffilter: \$(echo \$( vcffilter -h 2>&1) | grep 'vcflib' | sed 's/ filter.*\$//g' | sed 's/.* //g' )
                perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
            END_VERSIONS
            """
        }
        else {
            """
            echo tumor:$bams
            touch freebayes_${bed}.vcf

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                freebayes: \$(echo \$(freebayes --version 2>&1) | sed 's/version:\s*v//g' )
                vcffilter: \$(echo \$( vcffilter -h 2>&1) | grep 'vcflib' | sed 's/ filter.*\$//g' | sed 's/.* //g' )
                perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
            END_VERSIONS
            """
        }
}


/* Notes for freebayes :
freebayes \
  -t bedpart.9.bed \
  -f GCA_000001405.15_GRCh38_no_alt_analysis_set_nochr.fna \
  --pooled-continuous \
  --min-repeat-entropy 1 \
  -F 0.01 \
  -C 3 \
  -m 20 \
  -q 20 \
  --use-best-n-alleles 4 \
  GMSAL-AML5-T.T.sort.dedup.realign.bam \
  GMSAL-AML5-N.N.sort.dedup.realign.bam \
  > freebayes_bedpart.9.bed.vcf

*/