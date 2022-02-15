include { vep; snpeff } from '../process/annotation'

workflow ANNOTATION_TUMOR {
    take:
        vcfFile
        sampleName

    main:
        ANNOTATION_NORMAL(vcfFile, sampleName)

    emit:
        ANNOTATION_NORMAL.out
}

workflow ANNOTATION_NORMAL {
    take:
        vcfFile
        sampleName

    main:
        vep(vcfFile, sampleName)
        snpeff(vcfFile, sampleName)

    emit:
        vep.out.concat(snpeff.out)
}

