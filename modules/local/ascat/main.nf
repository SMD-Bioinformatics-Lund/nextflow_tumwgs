// ASCAT v3 allele-specific copy number for tumor/normal pairs (paired only).
// Split in two so a refit (ext.args of ASCAT_FIT) reuses the allele counts on -resume.
// Input tuples are grouped per case: meta/cram/crai/bai are lists, ordered by meta.type.

def ascat_sex(meta, idx) {
    def s = "${meta.sex[idx]}".toUpperCase()
    return (s in ['F', 'FEMALE', 'XX']) ? 'XX' : 'XY'
}

process ASCAT_PREPARE {
    label 'process_high'
    tag "${meta.id[tumor_idx]}"

    input:
        tuple val(group), val(meta), file(cram), file(crai), file(bai)

    output:
        tuple val(group), val(meta), path("*_{tumour,normal}{LogR,BAF}.txt"),    emit: logr_baf
        path "versions.yml",                                                    emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args   = task.ext.args ?: ''
        tumor_idx  = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
        normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N' }
        prefix = task.ext.prefix ?: "${meta.id[tumor_idx]}"
        sex    = ascat_sex(meta, normal_idx)
        """
        ascat.R prepare \\
            --tumor=${prefix} \\
            --normal=${meta.id[normal_idx]} \\
            --tumor_cram=${cram[tumor_idx]} \\
            --normal_cram=${cram[normal_idx]} \\
            --sex=${sex} \\
            --threads=${task.cpus} \\
            $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            ascat: \$(Rscript -e 'cat(as.character(packageVersion("ASCAT")))')
            allelecounter: \$(alleleCounter --version 2>&1 | head -1)
        END_VERSIONS
        """

    stub:
        tumor_idx  = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
        prefix = task.ext.prefix ?: "${meta.id[tumor_idx]}"
        """
        touch ${prefix}_tumourLogR.txt ${prefix}_tumourBAF.txt ${prefix}_normalLogR.txt ${prefix}_normalBAF.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            ascat: stub
            allelecounter: stub
        END_VERSIONS
        """
}

process ASCAT_FIT {
    label 'process_medium'
    tag "${meta.id[tumor_idx]}"

    input:
        tuple val(group), val(meta), path(logr_baf)

    output:
        tuple val(group), val(meta), path("*.ascat.purity_ploidy.tsv"), path("*.ascat.segments.tsv"), path("*.ascat.segments_raw.tsv"),  emit: segments
        tuple val(group), path("*.ascat.ownfit.*"), optional: true,                                   emit: ownfit   // ASCAT's own fit, after a both-allele refit
        tuple val(group), path("*.ascat.metrics.tsv"),                                                                                  emit: metrics
        tuple val(group), path("*.png"),                                                              emit: plots
        tuple val(group), path("*.ascat.RData"),                                                      emit: rdata
        path "versions.yml",                                                                          emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args   = task.ext.args ?: ''
        tumor_idx  = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
        normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N' }
        prefix = task.ext.prefix ?: "${meta.id[tumor_idx]}"
        sex    = ascat_sex(meta, normal_idx)
        """
        ascat.R fit \\
            --tumor=${prefix} \\
            --sex=${sex} \\
            $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            ascat: \$(Rscript -e 'cat(as.character(packageVersion("ASCAT")))')
        END_VERSIONS
        """

    stub:
        tumor_idx  = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
        prefix = task.ext.prefix ?: "${meta.id[tumor_idx]}"
        """
        printf 'sample\\tpurity\\tploidy\\tgoodness_of_fit\\tnote\\n${prefix}\\t0.8\\t2.3\\t98\\tstub\\n' > ${prefix}.ascat.purity_ploidy.tsv
        printf 'sample\\tchr\\tstartpos\\tendpos\\tnMajor\\tnMinor\\n${prefix}\\t1\\t1\\t248000000\\t1\\t1\\n' > ${prefix}.ascat.segments.tsv
        printf 'sample\\tchr\\tstartpos\\tendpos\\tnMajor\\tnMinor\\tnAraw\\tnBraw\\n${prefix}\\t1\\t1\\t248000000\\t1\\t1\\t1.0\\t1.0\\n' > ${prefix}.ascat.segments_raw.tsv
        touch ${prefix}.ascat.metrics.tsv ${prefix}.ASCATprofile.png ${prefix}.ascat.RData

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            ascat: stub
        END_VERSIONS
        """
}

process ASCAT_LOH_HEH {
    label 'process_single'
    tag "${meta.id[tumor_idx]}"

    input:
        tuple val(group), val(meta), path(purity_ploidy), path(segments), path(segments_raw)

    output:
        tuple val(group), val(meta), path("*.ascat.loh.tsv"), path("*.ascat.subclonal.tsv"), path("*.ascat.gains.tsv"), path("*.ascat.chrom.tsv"), path("*.ascat.heh.txt"),  emit: loh_heh
        tuple val(group), val(meta), path("*.ascat.cn_summary.json"), path("*.ascat.arms.tsv"), path("*.ascat.events.tsv"),                      emit: cn_summary
        path "versions.yml",                                                                                      emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args   = task.ext.args ?: ''
        tumor_idx  = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
        normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N' }
        prefix = task.ext.prefix ?: "${meta.id[tumor_idx]}"
        sex    = ascat_sex(meta, normal_idx)
        """
        cn_summary.py \\
            --purity_ploidy ${purity_ploidy} \\
            --segments ${segments} \\
            --segments_raw ${segments_raw} \\
            --prefix ${prefix} \\
            --label ascat \\
            --sex ${sex} \\
            $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python3 --version 2>&1 | sed 's/Python //')
        END_VERSIONS
        """

    stub:
        tumor_idx  = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
        prefix = task.ext.prefix ?: "${meta.id[tumor_idx]}"
        """
        touch ${prefix}.ascat.loh.tsv ${prefix}.ascat.subclonal.tsv ${prefix}.ascat.gains.tsv ${prefix}.ascat.chrom.tsv ${prefix}.ascat.heh.txt
        touch ${prefix}.ascat.cn_summary.json ${prefix}.ascat.arms.tsv ${prefix}.ascat.events.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: stub
        END_VERSIONS
        """
}
