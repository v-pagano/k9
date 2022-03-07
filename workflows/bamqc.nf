include { pb_collectmetrics } from '../process/parabricks'

workflow BAMQC {
    
    take:
        bam

    main:
        pb_collectmetrics(bam, params.pb_reference)

    emit:
        pb_collectmetrics.out

}