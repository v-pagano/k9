include { vep; snpeff } from '../process/annotation'

workflow ANNOTATION {
    take:
        vcf
        
    main:

        publishFiles = Channel.empty()

        vep(vcf, params.reference)
        publishFiles = publishFiles.mix(vep.out)

        snpeff(vcf)
        publishFiles = publishFiles.mix(snpeff.out)

    emit:
        publishFiles
}

