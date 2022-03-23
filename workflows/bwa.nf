include { bwa_mem; bwa_mem2 } from '../process/bwa'
include { SAM_SORT; SAM_MERGE; SAM_INDEX } from './sam'

workflow BWA {
    take: 
        fastq

    main:
        tempFiles = bwa_mem(fastq)
        unMergedFiles = SAM_SORT(tempFiles)
        SAM_MERGE(unMergedFiles, fastq, 'bwa')

    emit:
        bam = SAM_MERGE.out.bam
        publishFiles = SAM_MERGE.out.publishFiles
}

workflow BWA2 {
    take: 
        fastq

    main:
        tempFiles = bwa_mem2(fastq)
        unMergedFiles = SAM_SORT(tempFiles)
        SAM_MERGE(unMergedFiles, fastq, 'bwa2')

    emit:
        bam = SAM_MERGE.out.bam
        publishFiles = SAM_MERGE.out.publishFiles

}
