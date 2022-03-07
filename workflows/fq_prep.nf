include { fastqc } from '../process/fq_prep'

workflow FQ_PREP {

    take:
        fastq

    main:
        publishFiles = Channel.empty()

        fastqc(fastq)
        publishFiles = publishFiles.mix(fastqc.out[0], fastqc.out[1])

    emit:
        publishFiles   
}