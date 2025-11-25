process COYOTE {
    label "process_single"
    tag "$group"

    input:
        tuple val(group), val(meta), file(vcf), file(importy)

    output:
        tuple val(group), file("*.coyote"),  emit: coyote_import

    when:
        task.ext.when == null || task.ext.when

    script:

        process_group = group
        tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
        normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N' }
        if( meta.id.size() < 2 ) {
            process_group = group + '_single'
            tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
        }
        cnvplot_idx    = importy.findIndexOf{ it =~ 'modeled.png' }
        cnvseg_idx     = importy.findIndexOf{ it =~ 'panel' }
        fusions_idx    = importy.findIndexOf{ it =~ 'fusions' }

        fusions  = importy[fusions_idx]
        tumPlot  = importy[cnvplot_idx]
        cnv      = importy[cnvseg_idx]
        coyote_Group = params.coyote_group.split('-')[0]
        sample_subpanel = params.coyote_group.split('-')[1]

        """
        echo "/data/bnf/dev/saile/other/prj/import_scripts/import_myeloid_to_coyote_vep_gms_dev_WGS.pl \\
            --group ${coyote_Group} \\
            --subpanel ${sample_subpanel} \\
            --id ${process_group} \\
            --clarity-sample-id ${meta.clarity_sample_id[tumor_idx]} \\
            --clarity-pool-id ${meta.clarity_pool_id[tumor_idx]} \\
            --build 38 \\
            --vcf /access/${params.subdir}/vcf/${vcf} \\
            --cnv /access/tumwgs/cnv/${cnv} \\
            --cnvprofile  /access/tumwgs/cov/${tumPlot} \\
            --transloc /access/tumwgs/vcf/${fusions} \\
            --gens ${meta.id[tumor_idx]} \\
            --gensNorm ${meta.id[normal_idx]} " > ${process_group}.coyote
        """

    stub:
        process_group = group
        tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
        normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N' }
        if( meta.id.size() < 2 ) {
            process_group = group + '_single'
            tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
        }

        cnvplot_idx    = importy.findIndexOf{ it =~ 'modeled.png' }
        cnvseg_idx     = importy.findIndexOf{ it =~ 'panel' }
        fusions_idx    = importy.findIndexOf{ it =~ 'fusions' }

        fusions  = importy[fusions_idx]
        tumPlot  = importy[cnvplot_idx]
        cnv      = importy[cnvseg_idx]
        coyote_Group = params.coyote_group.split('-')[0]
        sample_subpanel = params.coyote_group.split('-')[1]
        """
        echo "/data/bnf/dev/saile/other/prj/import_scripts/import_myeloid_to_coyote_vep_gms_dev_WGS.pl \\
            --group ${coyote_Group} \\
            --subpanel ${sample_subpanel} \\
            --id ${process_group} \\
            --clarity-sample-id ${meta.clarity_sample_id[tumor_idx]} \\
            --clarity-pool-id ${meta.clarity_pool_id[tumor_idx]} \\
            --build 38 \\
            --vcf /access/${params.subdir}/vcf/${vcf} \\
            --cnv /access/tumwgs/cnv/${cnv} \\
            --cnvprofile  /access/tumwgs/cov/${tumPlot} \\
            --transloc /access/tumwgs/vcf/${fusions} \\
            --gens ${meta.id[tumor_idx]} \\
            --gensNorm ${meta.id[normal_idx]} " > ${process_group}.coyote
        """
}

process COYOTE_YAML {
    label "process_single"
    tag "$group"

    input:
       tuple val(group), val(meta), file(vcf), file(importy)

    output:
        tuple val(group), file("*.coyote3.yaml"), emit: coyote_import

    when:
        task.ext.when == null || task.ext.when

    script:
        environment = params.dev ? 'development' : params.validation ? 'validation' : params.testing ? 'testing' : 'production'
        process_group = group
        tumor_idx = 0
                tumor_sample = meta.id[tumor_idx]
        tumor_reads = meta.reads[tumor_idx] ?: null
        tumor_ffpe = meta.ffpe[tumor_idx] ? true : false
        tumor_sequencing_run = meta.sequencing_run[tumor_idx] ?: null
        tumor_purity = meta.purity[tumor_idx] ? meta.purity[tumor_idx].toFloat() : null
        normal_sample = null
        clarity_control_id = null
        clarity_control_pool_id = null
        control_reads = null
        control_ffpe = null
        control_sequencing_run = null
        control_purity = null
        sample_no = meta.id.size()
        coyote_Group = params.coyote_group.split('-')[0]
        sample_subpanel = params.coyote_group.split('-')[1]

        if( meta.id.size() >= 2 ) {
            process_group = group + 'p'
            tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
            normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N' }
            normal_sample = meta.id[normal_idx]
            sample_no = meta.id.size()
            clarity_control_id = meta.clarity_sample_id[normal_idx]
            clarity_control_pool_id = meta.clarity_pool_id[normal_idx]
            control_reads = meta.reads[normal_idx] ?: null
            control_ffpe = meta.ffpe[normal_idx] ? true : false
            control_sequencing_run = meta.sequencing_run[normal_idx] ?: null
            control_purity = meta.purity[normal_idx] ? meta.purity[normal_idx].toFloat() : null
        }

        cnvplot_idx    = importy.findIndexOf{ it =~ 'modeled.png' }
        cnvseg_idx     = importy.findIndexOf{ it =~ 'panel' }
        fusions_idx    = importy.findIndexOf{ it =~ 'fusions' }
        fusions  = importy[fusions_idx]
        tumPlot  = importy[cnvplot_idx]
        cnv      = importy[cnvseg_idx]

        """
        echo --- > ${process_group}.coyote3.yaml
        echo name: \\'${process_group}\\' >> ${process_group}.coyote3.yaml
        echo assay: \\'$params.coyote_group\\' >> ${process_group}.coyote3.yaml
        echo subpanel: \\'${meta.diagnosis[tumor_idx]}\\' >> ${process_group}.coyote3.yaml 
        echo sequencing_scope: \\'WGS\\' >> ${process_group}.coyote3.yaml
        echo omics_layer: \\'DNA\\' >> ${process_group}.coyote3.yaml
        echo sequencing_technology: \\'Illumina\\' >> ${process_group}.coyote3.yaml
        echo paired: ${meta.id.size() >= 2} >> ${process_group}.coyote3.yaml
        echo profile: \\'${environment}\\' >> ${process_group}.coyote3.yaml
        echo sample_no: ${sample_no} >> ${process_group}.coyote3.yaml
        echo case_id: \\'${tumor_sample}\\' >> ${process_group}.coyote3.yaml
        echo case_ffpe: ${tumor_ffpe} >> ${process_group}.coyote3.yaml
        echo case_sequencing_run: \\'${tumor_sequencing_run}\\' >> ${process_group}.coyote3.yaml
        echo case_reads: ${tumor_reads} >> ${process_group}.coyote3.yaml
        echo case_purity: ${tumor_purity} >> ${process_group}.coyote3.yaml
        echo control_id: \\'${normal_sample}\\' >> ${process_group}.coyote3.yaml
        echo clarity_case_pool_id: \\'${meta.clarity_pool_id[tumor_idx]}\\' >> ${process_group}.coyote3.yaml
        echo clarity_control_pool_id: \\'${clarity_control_pool_id}\\' >> ${process_group}.coyote3.yaml
        echo control_sequencing_run: \\'${control_sequencing_run}\\' >> ${process_group}.coyote3.yaml
        echo control_reads: ${control_reads} >> ${process_group}.coyote3.yaml
        echo control_purity: ${control_purity} >> ${process_group}.coyote3.yaml
        echo control_ffpe: ${control_ffpe} >> ${process_group}.coyote3.yaml
        echo genome_build: 38 >> ${process_group}.coyote3.yaml
        echo vcf_files: /access/${params.subdir}/vcf/${vcf} >> ${process_group}.coyote3.yaml
        echo cnvprofile: /access/${params.subdir}/cov/${tumPlot} >> ${process_group}.coyote3.yaml
        echo transloc: /access/${params.subdir}/manta/${fusions} >> ${process_group}.coyote3.yaml
        echo cnv: /access/${params.subdir}/cnv/${cnv} >> ${process_group}.coyote3.yaml
        echo pipeline: \\'${workflow.manifest.name}\\' >> ${process_group}.coyote3.yaml
        echo pipeline_version: ${workflow.manifest.version} >> ${process_group}.coyote3.yaml
        """
    stub:
        environment = params.dev ? 'development' : params.validation ? 'validation' : params.testing ? 'testing' : 'production'
        process_group = group
        tumor_idx = 0
                tumor_sample = meta.id[tumor_idx]
        tumor_reads = meta.reads[tumor_idx] ?: null
        tumor_ffpe = meta.ffpe[tumor_idx] ? true : false
        tumor_sequencing_run = meta.sequencing_run[tumor_idx] ?: null
        tumor_purity = meta.purity[tumor_idx] ? meta.purity[tumor_idx].toFloat() : null
        normal_sample = null
        clarity_control_id = null
        clarity_control_pool_id = null
        control_reads = null
        control_ffpe = null
        control_sequencing_run = null
        control_purity = null
        sample_no = meta.id.size()
        coyote_Group = params.coyote_group.split('-')[0]
        sample_subpanel = params.coyote_group.split('-')[1]

        if( meta.id.size() >= 2 ) {
            process_group = group + 'p'
            tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
            normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N' }
            normal_sample = meta.id[normal_idx]
            sample_no = meta.id.size()
            clarity_control_id = meta.clarity_sample_id[normal_idx]
            clarity_control_pool_id = meta.clarity_pool_id[normal_idx]
            control_reads = meta.reads[normal_idx] ?: null
            control_ffpe = meta.ffpe[normal_idx] ? true : false
            control_sequencing_run = meta.sequencing_run[normal_idx] ?: null
            control_purity = meta.purity[normal_idx] ? meta.purity[normal_idx].toFloat() : null
        }

        cnvplot_idx    = importy.findIndexOf{ it =~ 'modeled.png' }
        cnvseg_idx     = importy.findIndexOf{ it =~ 'panel' }
        fusions_idx    = importy.findIndexOf{ it =~ 'fusions' }
        fusions  = importy[fusions_idx]
        tumPlot  = importy[cnvplot_idx]
        cnv      = importy[cnvseg_idx]

        """
        echo --- > ${process_group}.coyote3.yaml
        echo name: \\'${process_group}\\' >> ${process_group}.coyote3.yaml
        echo assay: \\'$params.coyote_group\\' >> ${process_group}.coyote3.yaml
        echo subpanel: \\'${meta.diagnosis[tumor_idx]}\\' >> ${process_group}.coyote3.yaml 
        echo sequencing_scope: \\'WGS\\' >> ${process_group}.coyote3.yaml
        echo omics_layer: \\'DNA\\' >> ${process_group}.coyote3.yaml
        echo sequencing_technology: \\'Illumina\\' >> ${process_group}.coyote3.yaml
        echo paired: ${meta.id.size() >= 2} >> ${process_group}.coyote3.yaml
        echo profile: \\'${environment}\\' >> ${process_group}.coyote3.yaml
        echo sample_no: ${sample_no} >> ${process_group}.coyote3.yaml
        echo case_id: \\'${tumor_sample}\\' >> ${process_group}.coyote3.yaml
        echo case_ffpe: ${tumor_ffpe} >> ${process_group}.coyote3.yaml
        echo case_sequencing_run: \\'${tumor_sequencing_run}\\' >> ${process_group}.coyote3.yaml
        echo case_reads: ${tumor_reads} >> ${process_group}.coyote3.yaml
        echo case_purity: ${tumor_purity} >> ${process_group}.coyote3.yaml
        echo control_id: \\'${normal_sample}\\' >> ${process_group}.coyote3.yaml
        echo clarity_case_pool_id: \\'${meta.clarity_pool_id[tumor_idx]}\\' >> ${process_group}.coyote3.yaml
        echo clarity_control_pool_id: \\'${clarity_control_pool_id}\\' >> ${process_group}.coyote3.yaml
        echo control_sequencing_run: \\'${control_sequencing_run}\\' >> ${process_group}.coyote3.yaml
        echo control_reads: ${control_reads} >> ${process_group}.coyote3.yaml
        echo control_purity: ${control_purity} >> ${process_group}.coyote3.yaml
        echo control_ffpe: ${control_ffpe} >> ${process_group}.coyote3.yaml
        echo genome_build: 38 >> ${process_group}.coyote3.yaml
        echo vcf_files: /access/${params.subdir}/vcf/${vcf} >> ${process_group}.coyote3.yaml
        echo cnvprofile: /access/${params.subdir}/cov/${tumPlot} >> ${process_group}.coyote3.yaml
        echo transloc: /access/${params.subdir}/manta/${fusions} >> ${process_group}.coyote3.yaml
        echo cnv: /access/${params.subdir}/cnv/${cnv} >> ${process_group}.coyote3.yaml
        echo pipeline: \\'${workflow.manifest.name}\\' >> ${process_group}.coyote3.yaml
        echo pipeline_version: ${workflow.manifest.version} >> ${process_group}.coyote3.yaml
        """
}