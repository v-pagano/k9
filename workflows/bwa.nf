include { bwa_mem; bwa_mem2 } from '../process/bwa'
include { SAM_SORT; SAM_MERGE; SAM_INDEX } from './sam'

workflow BWA {
    take: 
        fastq

    main:
        tempFiles = bwa_mem(fastq)
        unMergedFiles = SAM_SORT(tempFiles)
        mergedFile = SAM_MERGE(unMergedFiles, fastq)

    emit:
        mergedFile
}

workflow BWA2 {
    take: 
        fastq

    main:
        tempFiles = bwa_mem2(fastq)
        unMergedFiles = SAM_SORT(tempFiles)
        mergedFile = SAM_MERGE(unMergedFiles, fastq)

    emit:
        mergedFile

}
