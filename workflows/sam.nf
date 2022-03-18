include { sam_sort; sam_merge } from '../process/sam'

workflow SAM_SORT {
    take: fileList

    main:
        sam_sort(fileList)

    emit:
        sam_sort.out
}

workflow SAM_MERGE {
    take: 
        fileList
        fastq

    main:
        sam_merge(fileList.toList(), fastq)

    emit:
        bam = sam_merge.out.bam
        publishFiles = sam_merge.out.publishFiles
}

workflow SAM_INDEX {
    take: 
        bamFile
        sampleName

    main:
        sam_index(bamFile, sampleName)

    emit:
        sam_index.out
}
