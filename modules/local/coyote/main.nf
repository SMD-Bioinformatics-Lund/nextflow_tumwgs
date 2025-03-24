process COYOTE {
    label "process_single"
    tag "$group"

    input:
        tuple val(group), val(meta), file(vcf), file(importy), file (cnvjson)

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

        def fusions = importy[0]
        def tumPlot = importy[1]
        def cnv = importy[2]
        def coyote_Group = params.coyote_group.split('-')[0]

        """
        echo "/data/bnf/scripts/import_myeloid_to_coyote_vep_gms_dev_WGS.pl --group ${coyote_Group} --id ${process_group} --vcf /access/${params.subdir}/vcf/${vcf} --cnv /access/tumwgs/cnv/${cnv} --transloc /access/tumwgs/manta/${fusions} --clarity-sample-id ${meta.clarity_sample_id[tumor_idx]} --build 38 --clarity-pool-id ${meta.clarity_pool_id[tumor_idx]} --gens ${meta.id[tumor_idx]} --cnvprofile /access/tumwgs/cov/${tumPlot}" > ${process_group}.coyote
        """

    stub:
        process_group = group
        tumor_idx = 0
        if( meta.id.size() < 2 ) {
            process_group = group + '_single'
            tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
        }

        def fusions = importy[0]
        def tumPlot = importy[1]
        def cnv = importy[2]
        def coyote_Group = params.coyote_group.split('-')[0]
        """
        
        echo "/data/bnf/scripts/import_myeloid_to_coyote_vep_gms_dev_WGS.pl --group ${coyote_Group} --id ${process_group} --vcf /access/${params.subdir}/vcf/${vcf} --cnv /access/${params.subdir}/cnv/${cnv} --transloc /access/${params.subdir}/manta/${fusions} --clarity-sample-id ${meta.clarity_sample_id[tumor_idx]} --build 38 --clarity-pool-id ${meta.clarity_pool_id[tumor_idx]} --gens ${meta.id[tumor_idx]} --cnvprofile /access/tumwgs/cov/${tumPlot}" > ${process_group}.coyote
        """
}

process COYOTE_YAML {
    label "process_single"
    tag "$group"

    input:
        tuple val(group), val(meta), file(vcf), file(importy),file (cnvjson)

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
        def fusions = importy[0]
        def tumPlot = importy[1]
        def cnv     = cnvjson[0]
        
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
        if( meta.id.size() < 2 ) {
            process_group = group + '_single'
            tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
        }
        def fusions = importy[0]
        def tumPlot = importy[1]
        def cnv     = cnvjson[0]
        
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