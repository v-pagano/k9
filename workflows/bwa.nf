include { bwa_mem; bwa_mem2 } from '../process/bwa'
include { SAM_SORT; SAM_MERGE; SAM_INDEX } from './sam'

workflow BWA {
    take: 
        fileTuple
        sampleName

    main:
        tempFiles = bwa_mem(fileTuple)
        unMergedFiles = SAM_SORT(tempFiles)
        mergedFile = SAM_MERGE(unMergedFiles, sampleName)

    emit:
        mergedFile
}

workflow BWA2 {
    take: 
        fileTuple
        sampleName

    main:
        tempFiles = bwa_mem2(fileTuple)
        unMergedFiles = SAM_SORT(tempFiles)
        mergedFile = SAM_MERGE(unMergedFiles, sampleName)

    emit:
        mergedFile

}
