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
        tumor_idx = 0
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
        echo "/data/bnf/dev/saile/other/prj/import_scripts/import_myeloid_to_coyote_vep_gms_dev_WGS.pl --group ${coyote_Group} --subpanel ${sample_subpanel} --id ${process_group} --vcf /access/${params.subdir}/vcf/${vcf} --cnv /access/tumwgs/cnv/${cnv} --transloc /access/tumwgs/manta/${fusions}  --cnvprofile  /access/tumwgs/cov/${tumPlot} --clarity-sample-id ${meta.clarity_sample_id[tumor_idx]} --build 38 --clarity-pool-id ${meta.clarity_pool_id[tumor_idx]} --gens ${meta.id[tumor_idx]} --cnvprofile /access/tumwgs/cov/${tumPlot}" > ${process_group}.coyote
        """

    stub:
        process_group = group
        tumor_idx = 0
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
        
        echo "/data/bnf/dev/saile/other/prj/import_scripts/import_myeloid_to_coyote_vep_gms_dev_WGS.pl --group ${coyote_Group}  --subpanel ${sample_subpanel} --id ${process_group} --vcf /access/${params.subdir}/vcf/${vcf} --cnv /access/${params.subdir}/cnv/${cnv} --transloc /access/${params.subdir}/manta/${fusions} --clarity-sample-id ${meta.clarity_sample_id[tumor_idx]} --build 38 --clarity-pool-id ${meta.clarity_pool_id[tumor_idx]} --gens ${meta.id[tumor_idx]} --cnvprofile /access/tumwgs/cov/${tumPlot}" > ${process_group}.coyote
        """
}

process COYOTE_YAML {
    label "process_single"
    tag "$group"

    input:
       tuple val(group), val(meta), file(vcf), file(importy)

    output:
        tuple val(group), file("*.coyote.yaml"), emit: coyote_import

    when:
        task.ext.when == null || task.ext.when

    script:
        process_group = group
        tumor_idx = 0
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
        echo --- > ${process_group}.coyote.yaml
        echo groups: [\\'$params.coyote_group\\'] >> ${process_group}.coyote.yaml
        echo subpanel: \\'${meta.diagnosis[tumor_idx]}\\' >> ${process_group}.coyote.yaml
        echo name: \\'${process_group}\\' >> ${process_group}.coyote.yaml
        echo clarity-sample-id: \\'${meta.clarity_sample_id[tumor_idx]}\\' >> ${process_group}.coyote.yaml
        echo clarity-pool-id: \\'${meta.clarity_pool_id[tumor_idx]}\\' >> ${process_group}.coyote.yaml
        echo genome_build: 38 >> ${process_group}.coyote.yaml
        echo vcf_files: /access/${params.subdir}/vcf/${vcf} >> ${process_group}.coyote.yaml
        echo cnvprofile: /access/${params.subdir}/cov/${tumPlot} >> ${process_group}.coyote.yaml
        echo transloc: /access/${params.subdir}/manta/${fusions} >> ${process_group}.coyote.yaml
        echo cnv: /access/${params.subdir}/cnv/${cnv} >> ${process_group}.coyote.yaml
        """
    stub:
        process_group = group
        tumor_idx = 0
        if( meta.id.size() <, 2 ) {
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
        echo --- > ${process_group}.coyote.yaml
        echo groups: [\\'$params.coyote_group\\'] >> ${process_group}.coyote.yaml
        echo subpanel: \\'${meta.diagnosis[tumor_idx]}\\' >> ${process_group}.coyote.yaml
        echo name: \\'${process_group}\\' >> ${process_group}.coyote.yaml
        echo clarity-sample-id: \\'${meta.clarity_sample_id[tumor_idx]}\\' >> ${process_group}.coyote.yaml
        echo clarity-pool-id: \\'${meta.clarity_pool_id[tumor_idx]}\\' >> ${process_group}.coyote.yaml
        echo genome_build: 38 >> ${process_group}.coyote.yaml
        echo vcf_files: /access/${params.subdir}/vcf/${vcf} >> ${process_group}.coyote.yaml
        echo cnvprofile: /access/${params.subdir}/cov/${tumPlot} >> ${process_group}.coyote.yaml
        echo transloc: /access/${params.subdir}/manta/${fusions} >> ${process_group}.coyote.yaml
        echo cnv: /access/${params.subdir}/cnv/${cnv} >> ${process_group}.coyote.yaml
        """
}