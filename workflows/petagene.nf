include { bam2pgbam } from '../process/petagene'

workflow BAM2PGBAM {
    take: 
        bam

    main:
        bam2pgbam(bam)

    emit:
        bam2pgbam.out
}