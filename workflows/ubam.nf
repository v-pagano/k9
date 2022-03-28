include { sam_bam2fq } from '../process/sam'
include { picard_fq2ubam } from '../process/picard'

workflow UBAM2FQ {
    take:
        bam

    main:
        sam_bam2fq(bam)

    emit:
        sam_bam2fq.out

}

workflow FQ2UBAM {
    take:
        fastq

    main:
        picard_fq2ubam(fastq)

    emit:
        picard_fq2ubam.out[1]

}