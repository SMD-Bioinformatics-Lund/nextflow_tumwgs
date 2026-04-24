#!/usr/bin/env nextflow

include { COYOTE               } from '../../modules/local/coyote/main'
include { COYOTE_YAML          } from '../../modules/local/coyote/main'

workflow ADD_TO_DB {
    take: 
        vcf             // channel: [mandatory] [ val(group), val(meta), file(vcf) ]
        cnvbed          // channel: [optional] [ val(group), file(cnvbed) ] 
        cnvjson            // channel: [optional] [ val(group), file(cnvjson) ]           
        fusions         // channel: [optional] [ val(group), file(segments) ] // have to change this
        tum_plot        // channel: [optional] [ val(group), file(segments) ]

    main:
        optional = cnvbed.mix(fusions,tum_plot).groupTuple()
        optional_json = cnvjson.mix(fusions,tum_plot).groupTuple()

        // optional.view()
        // vcf.view()
        // vcf.join(optional).view()
        COYOTE { vcf.join(optional) }
        COYOTE_YAML { vcf.join(optional_json) }

    emit:
        coyotedone = COYOTE.out.coyote_import        // channel: [ val(group), file(coyote) ]
        
}
