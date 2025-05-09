process PON_FILTER {
    label "process_medium"
    tag "$group"

    input:
        tuple val(group), val(meta), file(vcf) 
        
    output:
        tuple val(group), val(meta), file("*.agg.pon.vcf"), emit: vcf_pon
        path "versions.yml",                                emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix = task.ext.prefix ?: "${group}"
        def pons = []
        if( params.freebayes ) { pons.push("freebayes="+params.PON_freebayes) }
        if( params.vardict )   { pons.push("vardict="+params.PON_vardict) }
        //if( params.tnscope )   { pons.push("tnscope="+params.PON_tnscope) }
        def pons_str = pons.join(",")
        tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
        """
        filter_with_pon.pl --vcf $vcf --pons $pons_str --tumor-id ${meta.id[tumor_idx]} > ${prefix}.agg.pon.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${group}"
        tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
        """
        echo ${meta.id[tumor_idx]}
        touch ${prefix}.agg.pon.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
        END_VERSIONS
        """
}

process ANNOTATE_VEP {

    label "process_medium"
    tag "$group"

    input:
        tuple val(group), val(meta), file(vcf)
        
    output:
        tuple val(group), val(meta), file("*.vep.vcf"), emit: vcf_vep
        path "versions.yml",                            emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args   = task.ext.args   ?: ''
        def prefix = task.ext.prefix ?: "${vcf.baseName}"

        """
        vep -i ${vcf} -o ${prefix}.vep.vcf --fork ${task.cpus} $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            vep: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : //;s/ .*\$//')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${vcf.baseName}"

        """
        touch ${prefix}.vep.vcf
        echo $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
        END_VERSIONS
        """ 
}

process FILTER_PANEL{

    label "process_single"
    tag "$group"

    input:
        tuple val(group), val(meta), file(vcf) 

    output:
        tuple val(group), val(meta), file("*.agg.pon.vep.panel.vcf"),   emit: vcf_panel
        path "versions.yml",                                            emit: versions

    
    when:
        task.ext.when == null || task.ext.when

    script:

        def should_hard_filter = params.SNV_HARD_FILTER ? '1' : ''
        def args   = task.ext.args   ?: ''
        def prefix = task.ext.prefix ?: "${vcf.baseName}"

        """
        filter_with_panel_snv.pl $vcf $args $should_hard_filter > ${prefix}.panel.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${group}"
        tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
        """
        echo ${meta.id[tumor_idx]}
        touch ${prefix}.agg.pon.vep.panel.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
        END_VERSIONS
        """
}

process MARK_GERMLINES {
    label "process_single"
    tag "$group"

    input:
        tuple val(group), val(meta), file(vcf) // from vcf_germline.join(meta_germline.groupTuple())

        
    output:
        tuple val(group), val(meta), file("*.agg.pon.vep.markgerm.vcf"),    emit: vcf_germline
        path "versions.yml",                                                emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: "${group}"
        if( meta.id.size() >= 2 ) {
            tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
            normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N' }
            """
            fix_vep_gnomad.pl $vcf > ${prefix}.agg.pon.vep.fix.vcf
            mark_germlines.pl --vcf ${prefix}.agg.pon.vep.fix.vcf --tumor-id ${meta.id[tumor_idx]} --normal-id ${meta.id[normal_idx]} $args > ${prefix}p.agg.pon.vep.markgerm.vcf

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
            END_VERSIONS
            """
        }
        else if( meta.id.size() == 1 ) {
            """
            fix_vep_gnomad.pl $vcf > ${prefix}.agg.pon.vep.fix.vcf
            mark_germlines.pl --vcf ${prefix}.agg.pon.vep.fix.vcf --tumor-id ${meta.id[0]} $args > ${prefix}.agg.pon.vep.markgerm.vcf

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
            END_VERSIONS
            """
        }

    stub:
        def prefix = task.ext.prefix ?: "${group}"
        if( meta.id.size() >= 2 ) {
            tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
            normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N' }
            """
            echo --tumor-id ${meta.id[tumor_idx]} --normal-id ${meta.id[normal_idx]}
            touch ${prefix}p.agg.pon.vep.markgerm.vcf

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
            END_VERSIONS
            """
        }
        else if( meta.id.size() == 1 ) {
            """
            echo ${meta.id[0]}
            touch ${prefix}.agg.pon.vep.markgerm.vcf

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
            END_VERSIONS
            """
        }
}

process FILTER_FOR_CNV {
    label "process_single"
    tag "$group"

    input:
        tuple val(group), val(meta), file(vcf), val(vc), file(vcf_unfilt)

    output:
        tuple val(group), file("*_vardict.germlines.vcf.gz"), file("*_vardict.germlines.vcf.gz.tbi"),   emit: vcf_only_germline
        path "versions.yml",                                                                            emit: versions
    
    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: "${group}"
        """
        germline_for_cnvkit.pl $vcf > ${prefix}.agg.pon.vep.germline.vcf
        bedtools intersect -a $vcf_unfilt -b ${prefix}.agg.pon.vep.germline.vcf $args > ${prefix}_vardict.germlines.vcf
        bgzip ${prefix}_vardict.germlines.vcf
        tabix ${prefix}_vardict.germlines.vcf.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
            bgzip: \$(bgzip --v | grep 'bgzip' | sed 's/.* //g')
            tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
            perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${group}"
        """
        echo $vcf $vcf_unfilt
        touch ${prefix}_vardict.germlines.vcf.gz 
        touch ${prefix}_vardict.germlines.vcf.gz.tbi

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
            bgzip: \$(bgzip --v | grep 'bgzip' | sed 's/.* //g')
            tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
            perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
        END_VERSIONS
        """
}

process COYOTE_SEGMENTS {
    label "process_single"
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), file(vcf)
    
    output:
        tuple val(group), val(meta), file("*.cn-segments.panel.bed"),  emit: filtered
        tuple val(group), val(meta), file("*.cn-segments.bed"),        emit: raw
        path "versions.yml",                                           emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args   = task.ext.args   ?: ''
        def prefix = task.ext.prefix ?: "${meta.id}"
        def normal = meta.type.equals('normal') || meta.type.equals('N') ? "--normal" : ""

        """
        coyote_segmentator.pl --vcf $vcf --id ${meta.id} $normal $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        touch ${prefix}.cn-segments.panel.bed ${prefix}.cn-segments.bed

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
        END_VERSIONS
        """
}

process COYOTE_SEGMENTS_JSON {
    label "process_single"
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), file(bed)
    
    output:
        tuple val(group), val(meta), file("*panelmatched.json"),  emit: json_panel
        path "versions.yml",                                      emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args   = task.ext.args   ?: ''
        def prefix = task.ext.prefix ?: "${meta.id}"
        def normal = meta.type.equals('normal') || meta.type.equals('N') ? "--normal" : ""

        """
        cnvJSON.py --bed $bed $args --id ${meta.id} $normal

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version 2>&1| sed -e 's/Python //g')
        END_VERSIONS
        """

    stub:
        def args   = task.ext.args   ?: ''
        def prefix = task.ext.prefix ?: "${meta.id}"
        def normal = meta.type.equals('normal') || meta.type.equals('N') ? "--normal" : ""
        """
        echo cnvJSON.py --bed $bed $args --id ${meta.id} $normal > ${meta.id}cnvs_panelmatched.json

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version 2>&1| sed -e 's/Python //g')
        END_VERSIONS
        """
}

process MERGE_SEGMENTS {
    label "process_single"
    tag "$group"

    input:
        tuple val(group), val(meta), file(segments)

    output:
        tuple val(group), file("*.cn-segments.panel.merged.bed"), emit: merged

    when:
        task.ext.when == null || task.ext.when

    script:
        def args   = task.ext.args   ?: ''
        def prefix = task.ext.prefix ?: "${group}"
        """
        cat $segments > ${prefix}.cn-segments.panel.merged.bed
        """

    stub:
        def prefix = task.ext.prefix ?: "${group}"
        """
        touch ${prefix}.cn-segments.panel.merged.bed
        """

}

process MERGE_JSON {
    label "process_single"
    tag "$group"

    input:
        tuple val(group), val(meta), file(segments)

    output:
        tuple val(group), file("*.cnvs.merged.json"), emit: merged

    when:
        task.ext.when == null || task.ext.when

    script:
        def args   = task.ext.args   ?: ''
        def prefix = task.ext.prefix ?: "${group}"
        segments = segments.join(' ')
        if (meta.id.size() >= 2) {
            """
            jq $args $segments > ${group}.cnvs.merged.json
            """
        }
        else {
            """
            cat $segments > ${group}.cnvs.merged.json
            """
        }


    stub:
        def args   = task.ext.args   ?: ''
        def prefix = task.ext.prefix ?: "${group}"
        segments = segments.join(' ')
        if (meta.id.size() >= 2) {
            """
            echo jq $args $segments > ${group}.cnvs.merged.json
            """
        }
        else {
            """
            echo $segments > ${group}.cnvs.merged.json
            """
        }

}

process FILTER_MANTA {
    label "process_single"
    tag "$group"

    input:
        tuple val(group), val(meta), file(vcf)

    output:
        tuple val(group), val(meta), file("*_manta_filtered.vcf"),     emit: filtered
        tuple val(group), val(meta), file("*_manta_bnd_filtered.vcf"), emit: bnd_filtered
        path "versions.yml",                                           emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args   = task.ext.args   ?: ''
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        filter_manta.pl --vcf $vcf $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        touch ${prefix}_manta_bnd_filtered.vcf
        touch ${prefix}_manta_filtered.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
        END_VERSIONS
        """

}

process GENEFUSE_JSON_TO_VCF {
    label "process_single"
    tag "$group"
    
    input:
        tuple val(group), val(meta), file(json)

    output:
        tuple val(group), file("*_genefuse.vcf"),  emit: genefuse_vcf
        path "versions.yml",                                emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        genefuse_json_to_vcf.py -i ${meta.id} -j $json -o ${prefix}_genefuse.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version 2>&1| sed -e 's/Python //g')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        touch ${prefix}_genefuse.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version 2>&1 | sed -e 's/Python //g')
        END_VERSIONS
        """
}

process BIOMARKERS_TO_JSON {
    label "process_single"
    tag "$group"

    input:
        tuple val(group), file(markers)

    output:
        tuple val(group), file("*.bio.json"),   emit: biomarkers_json
        path "versions.yml",                    emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix = task.ext.prefix ?: "${group}"
        msis_idx = markers.findIndexOf{ it =~ 'msi_single' }
        msip_idx = markers.findIndexOf{ it =~ 'msi_paired' }
        hrd_idx = markers.findIndexOf{ it =~ 'HRD' }
        // find biomarkers //
        msis = msis_idx >= 0 ? markers[msis_idx].collect {'--msi_s ' + it} : null
        msip = msip_idx >= 0 ? markers[msip_idx].collect {'--msi_p ' + it} : null
        hrd = hrd_idx >= 0 ? markers[hrd_idx].collect {'--hrd ' + it} : null
        tmp = (msis ?: []) + (msip ?: []) + (hrd ?: [])
        command = tmp.join(' ')

        """
        aggregate_biomarkers.py $command --out ${prefix}.bio.json --id $group

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version 2>&1| sed -e 's/Python //g')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${group}"
        msis_idx = markers.findIndexOf{ it =~ 'msi_single' }
        msip_idx = markers.findIndexOf{ it =~ 'msi_paired' }
        hrd_idx = markers.findIndexOf{ it =~ 'HRD' }
        // find biomarkers //
        msis = msis_idx >= 0 ? markers[msis_idx].collect {'--msi_s ' + it} : null
        msip = msip_idx >= 0 ? markers[msip_idx].collect {'--msi_p ' + it} : null
        hrd = hrd_idx >= 0 ? markers[hrd_idx].collect {'--hrd ' + it} : null
        tmp = (msis ?: []) + (msip ?: []) + (hrd ?: [])
        command = tmp.join(' ')

        """
        echo $command
        touch ${prefix}.bio.json

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version 2>&1| sed -e 's/Python //g')
        END_VERSIONS
        """
}

process VCFANNO {
    label "process_single"
    tag "$group"

    input:
        tuple val(group), val(meta), file(vcf) 
        
    output:
        tuple val(group), val(meta), file("*.agg.enigma.vcf"),  emit: vcf_enigma
        path "versions.yml",                                    emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: "${group}"
        """
        vcfanno_linux64 $args $vcf > ${prefix}.agg.enigma.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            vcfanno: \$( echo \$(vcfanno_linux64 2>&1) | sed 's/.*version //' | sed 's/ \\[.*//')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${group}"
        """
        touch ${prefix}.agg.enigma.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            vcfanno: \$( echo \$(vcfanno_linux64 2>&1) | sed 's/.*version //' | sed 's/ \\[.*//')
        END_VERSIONS
        """
}

process CREATE_SNVPON {
    label "process_alot"
    tag "$vc"

    input:
        tuple val(group), val(vc), file(vcfs) 

    output:
        tuple val(group), val(vc), file("*_${vc}_PON.snv"), emit: SNV_PON
        path "versions.yml",                                emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix = task.ext.prefix ?: "${params.assay}"
        """
        create_snv_pon.pl "*.vcf.gz" > ${prefix}_${vc}_PON.snv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${params.assay}"
        """
        echo $vcfs
        touch ${prefix}_${vc}_PON.snv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
        END_VERSIONS
        """
}

process CONTAMINATION {
    label "process_single"
    tag "$group"

    input:
        tuple val(group), val(meta), file(vcf)

    output:
        tuple val(group), file("*.txt"), file("*.png"), emit: contamination_result_files
        tuple val(group), file("*.contaminationpy"),    emit: contamination_cdm
        path "versions.yml",                            emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args     ?: ''
        def args2   = task.ext.args2    ?: ''

        if(meta.id.size() >= 2) { 
            tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
            normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N' }

            """
            find_contaminant.pl --vcf $vcf --case-id ${meta.id[tumor_idx]} $args > ${meta.id[tumor_idx]}.value
            echo "--overwrite --sample-id ${meta.id[tumor_idx]} --run-folder ${meta.sequencing_run[tumor_idx]} --assay ${params.cdm} --contamination" > ${meta.id[tumor_idx]}.1
            paste -d " " ${meta.id[tumor_idx]}.1 ${meta.id[tumor_idx]}.value > ${meta.id[tumor_idx]}.contamination
            find_contaminant.pl --vcf $vcf --case-id ${meta.id[tumor_idx]} $args2 > ${meta.id[normal_idx]}.value
            echo "--overwrite --sample-id ${meta.id[normal_idx]} --sequencing-run ${meta.sequencing_run[normal_idx]} --assay ${params.cdm} --contamination" > ${meta.id[normal_idx]}.1
            paste -d " " ${meta.id[normal_idx]}.1 ${meta.id[normal_idx]}.value > ${meta.id[normal_idx]}.contaminationpy

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
            END_VERSIONS
            """
        }
        else {
            """
            find_contaminant.pl --vcf $vcf --case-id ${meta.id[0]} $args > ${meta.id[0]}.value
            echo "--overwrite --sample-id ${meta.id[0]} --sequencing-run ${meta.sequencing_run[0]} --assay ${params.cdm} --contamination" > ${meta.id[0]}.1
            paste -d " " ${meta.id[0]}.1 ${meta.id[0]}.value > ${meta.id[0]}.contaminationpy

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
            END_VERSIONS
            """
        }

    stub:
        if(meta.id.size() >= 2) { 
            tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
            normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N' }
            """
            touch test.png
            touch test.txt
            touch ${meta.id[tumor_idx]}.contaminationpy
            touch ${meta.id[normal_idx]}.contaminationpy

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
            END_VERSIONS
            """
        }
        else {
            """
            touch test.png
            touch test.txt
            touch ${meta.id[0]}.contaminationpy

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
            END_VERSIONS
            """
        }
}

process BEDTOOLS_INTERSECT {
    label "process_single"
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), val(vc), file(vcf)
        val(bed)

    output:
        tuple val(group), val(vc), file("*_${vc}_intersected.vcf"), emit: vcf_intersected
        path "versions.yml",                                        emit: versions
    
    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args     ?: ''
        def prefix  = task.ext.prefix ?: "${meta.id}"
        """
        bedtools intersect -a $vcf -b $bed $args > ${prefix}_${vc}_intersected.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
        END_VERSIONS
        """

    stub:
        def prefix  = task.ext.prefix ?: "${meta.id}"
        """
        echo $vcf
        touch ${prefix}_${vc}_intersected.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
        END_VERSIONS
        """
}

process OVERLAP_GENES {
    label "process_low"
    tag "${meta.group}"

    input:
		tuple val(group), val(meta), file(segments) 

    output:
        tuple val(group), val(meta), file("*.cnv.annotated.bed"),       emit : annotated_bed
        path "versions.yml",                                            emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T'  }
    def normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N'  }
    def prefix = task.ext.prefix ?: "${meta.id[tumor_idx]}"
        """
        overlapping_genes.pl ${segments} ${args} > ${prefix}.cnv.annotated.bed
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
        END_VERSIONS        
        """

    stub:
    def tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T'  }
    def normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N'  }
    def prefix = task.ext.prefix ?: "${meta.id[tumor_idx]}"
        """
        
        touch ${prefix}.cnv.annotated.bed
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
        END_VERSIONS        
        """
}

process FILTER_CNVS_PANEL {
    label "process_single"
    tag "$group"

    input:
        tuple val(group), val(meta), file(bed) 

    output:
        tuple val(group), file("*.cnv.annotated.panel.bed"), file("*.cnv.annotated.panel.json"),            emit: vcf_panel
        path "versions.yml",                                                                                emit: versions

    
    when:
        task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T'  }
    def normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N'  }
    def prefix = task.ext.prefix ?: "${meta.id[tumor_idx]}"
        """
        filter_with_panel_cnv.pl ${bed} ${args} > ${prefix}.cnv.annotated.panel.bed
        cnv_json.py --bed ${prefix}.cnv.annotated.panel.bed  --json ${prefix}.cnv.annotated.panel.json

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
        END_VERSIONS
        """

    stub:
    def tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T'  }
    def normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N'  }
    def prefix = task.ext.prefix ?: "${meta.id[tumor_idx]}"
        """
        touch ${prefix}.cnv.annotated.panel.bed
        touch ${prefix}.cnv.annotated.panel.json

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
        END_VERSIONS
        
        """
}

process FILTER_FUSIONS_PANEL {
    label "process_single"
    tag "$group"

    input:
        tuple val(group), val(meta), file(vcf) 

    output:
        tuple val(group), val(meta), file("*.manta.fusions.vcf"),               emit: sv_panel
        path "versions.yml",                                                    emit: versions

    
    when:
        task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${group}"
        """
        filter_with_panel_fusions.pl ${vcf} ${args} > ${prefix}.manta.fusions.vcf	

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
        END_VERSIONS
        """

    
    stub:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${group}"
        """
        echo $args
        touch ${prefix}.manta.fusions.vcf	

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
        END_VERSIONS
        """
}

process FIX_VEP {
    label "process_single"
    tag "$group"

    input:
        tuple val(group), val(meta), file(vcf) // from vcf_germline.join(meta_germline.groupTuple())

        
    output:
        tuple val(group), val(meta), file("*.agg.pon.vep.fix.vcf"),             emit: fixed_vcf
        path "versions.yml",                                                    emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: "${group}"
        if( meta.id.size() >= 2 ) {
            tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
            normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N' }
            """
            fix_vep_gnomad.pl $vcf > ${prefix}.agg.pon.vep.fix.vcf

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
            END_VERSIONS
            """
        }
        else if( meta.id.size() == 1 ) {
            """
            fix_vep_gnomad.pl $vcf > ${prefix}.agg.pon.vep.fix.vcf

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
            END_VERSIONS
            """
        }

    stub:
        def prefix = task.ext.prefix ?: "${group}"
        if( meta.id.size() >= 2 ) {
            tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
            normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N' }
            """
            echo --tumor-id ${meta.id[tumor_idx]} --normal-id ${meta.id[normal_idx]}
            touch ${prefix}p.agg.pon.vep.fix.vcf

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
            END_VERSIONS
            """
        }
        else if( meta.id.size() == 1 ) {
            """
            echo ${meta.id[0]}
            touch ${prefix}.agg.pon.vep.fix.vcf

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
            END_VERSIONS
            """
        }
}

process POST_ANNOTATION_FILTERS {
    label "process_single"
    tag "$group"

    input:
        tuple val(group), val(meta), file(vcf)
        
    output:
        tuple val(group), val(meta), file("*.final.filtered.vcf"),              emit: filtered_vcf
        path "versions.yml",                                                    emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix  = task.ext.prefix   ?: "${group}"
        def args    = task.ext.args     ?: ''
        """
        post_annotation_filtering.py --vcf $vcf $args > ${prefix}.final.filtered.vcf
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version 2>&1| sed -e 's/Python //g')
        END_VERSIONS
        """
    stub:
        def prefix  = task.ext.prefix ?: "${group}"
        def args    = task.ext.args     ?: ''
        """
        touch ${prefix}.final.filtered.vcf
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version 2>&1| sed -e 's/Python //g')
        END_VERSIONS
        """
        
}

