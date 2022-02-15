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
        sampleName

    main:
        sam_merge(fileList.toList(), sampleName)

    emit:
        sam_merge.out
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
