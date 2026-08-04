process CREATE_PED_FILES {
    label "process_low"
    tag "$group_id"

    input:
        tuple val(group_id), val (meta)

    output:
        tuple val(group_id), val(meta), path("${meta.group}_${meta.id}_base.ped"), emit: ped_file
        path "versions.yml",                                                       emit: versions
        
    script:
        def father = meta.father ?: "0"
        def mother = meta.mother ?: "0"
        def phenotype = meta.phenotype ?: "0"  

        """
        create_ped.pl --mother ${mother} --father ${father} --group ${meta.group} --id ${meta.id} --sex ${meta.sex}
        mv ${meta.group}_base.ped ${meta.group}_${meta.id}_base.ped

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
        END_VERSIONS
        """
    stub:
        def father = meta.father ?: "0"
        def mother = meta.mother ?: "0"
        def phenotype = meta.phenotype ?: "0"  

        """
        touch ${meta.group}_${meta.id}_base.ped

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
        END_VERSIONS
        """
}


process SOMALIER_QC {
    label "process_medium"
    tag "${group_id}"

    input:
        tuple val(group_id), val(meta), file(crams), file(crai), file(bai), file (ped_files)

    output:
        tuple val(group_id), val(meta), file("*.samples.tsv"), file("*.pairs.tsv"),  file("*.groups.tsv"), file("*.contamination.tsv"), emit: somalier_check
        path "versions.yml", emit: versions

    script:
        def args = task.ext.args ?: "" 
        def args2 = task.ext.args2 ?: ""  

        tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
        normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N' }
        tumor_bam = crams[tumor_idx]
        normal_bam = crams[normal_idx]
        tumor_id = meta.id[tumor_idx]
        normal_id = meta.id[normal_idx]
        ped_normal = ped_files[normal_idx]
        ped_tumor = ped_files[tumor_idx]

        
        """
        somalier extract -d extracted $args ${tumor_bam} 
        somalier extract -d extracted $args ${normal_bam}
        cat ${ped_normal} ${ped_tumor}  > ${group_id}.merged.ped
        somalier relate --ped ${group_id}.merged.ped --infer extracted/*.somalier $args2 -o ${group_id}   
        somalier contamination -p ./extracted/${tumor_id}.somalier ./extracted/${normal_id}.somalier $args2 -o ${group_id}.contamination.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            somalier: \$(somalier 2>&1 |sed -n 's/.*version: \\([0-9.]*\\).*/\\1/p')
        END_VERSIONS
        """

      stub:
        tumor_idx  = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
        normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N' }
        tumor_id   = meta.id[tumor_idx]
        normal_id  = meta.id[normal_idx]

        """
        touch ${group_id}.samples.tsv ${group_id}.pairs.tsv ${group_id}.contamination.tsv ${group_id}.groups.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            somalier: \$(somalier 2>&1 |sed -n 's/.*version: \\([0-9.]*\\).*/\\1/p')
        END_VERSIONS
        """
}


process SOMALIER2CDM {
    label "process_single"
    tag "${meta.id}"

    input: 
        tuple val(group), val(meta), file(samples_stats), file(pairs_stats),  file(groups), file(contamination)

    output:
        tuple val(group), val(meta), path("*somalier.json"), emit: json
        tuple val(group), val(meta), path("*peddy2cdm"), emit: cdm
        path "versions.yml", emit: versions

    script:
        def args    = task.ext.args     ?: ''
        def args2   = task.ext.args2    ?: ''

        def tumor_idx   = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
        def normal_idx  = meta.type.findIndexOf{ it == 'normal' || it == 'N' }
        def tumor_id = meta.id[tumor_idx]
        def normal_id = meta.id[normal_idx]
        def tumor_run = meta.sequencing_run[tumor_idx]
        def normal_run = meta.sequencing_run[normal_idx]

        def tumor_arg =  "${tumor_id}:${tumor_run}"
        def normal_arg =  "${normal_id}:${normal_run}"
            
        """
        somalier2json.py \
        --somalier $samples_stats \
        --sample $tumor_arg \
        --sample $normal_arg \
        $args $args2

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version 2>&1| sed -e 's/Python //g')
        END_VERSIONS
        """
        
    stub:
        def args    = task.ext.args     ?: ''
        def args2   = task.ext.args2    ?: ''

        def tumor_idx   = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
        def normal_idx  = meta.type.findIndexOf{ it == 'normal' || it == 'N' }
        def tumor_id = meta.id[tumor_idx]
        def normal_id = meta.id[normal_idx]
        def tumor_run = meta.sequencing_run[tumor_idx]
        def normal_run = meta.sequencing_run[normal_idx]

        def tumor_arg =  "${tumor_id}:${tumor_run}"
        def normal_arg =  "${normal_id}:${normal_run}"

        """
        echo "somalier2json.py --somalier $samples_stats --sample $tumor_arg --sample $normal_arg $args $args2"
        touch "${tumor_id}.somalier.json"
        touch "${tumor_id}.peddy2cdm"
        touch "${normal_id}.somalier.json"
        touch "${normal_id}.peddy2cdm"
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version 2>&1| sed -e 's/Python //g')
        END_VERSIONS
        """
}

