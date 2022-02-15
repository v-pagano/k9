include { pb_fq2bam; pb_germline; pb_haplotypecaller; pb_deepvariant; pb_somatic } from '../process/parabricks.nf'

workflow PB_FQ2BAM {
    take: 
        fileTuple
        sampleName

    main:
        infq = []
        infq = fileTuple.collect { ' --in-fq ' + it.fastq1 + ' ' + it.fastq2 + ' "' + it.rg + '" ' }

        pb_fq2bam(infq.join(' '), sampleName)

    emit:
        pb_fq2bam.out
}

workflow PB_RNA_FQ2BAM {
    take: fileTuple

    main:
        pb_rna_fq2bam(fileTuple.flatten().map { ' --in-fq ' + it.r1 + ' ' + it.r2 + ' "' + it.rg + '" ' }.toList())
    emit:
        pb_rna_fq2bam.out
}

workflow PB_SOMATIC {
    take:         
        fileNormal
        sampleNormal
        fileTumor
        sampleTumor

    main:
        nfq = fileNormal.collect { ' ' + it.fastq1 + ' ' + it.fastq2 + ' "' + it.rg + '" ' }
        tfq = fileTumor.collect { ' ' + it.fastq1 + ' ' + it.fastq2 + ' "' + it.rg + '" ' }
        pb_somatic(tfq.join(' '), nfq.join(' '), sampleTumor, sampleNormal)

    emit:
        pb_somatic.out
}

workflow PB_GERMLINE {
    take:         
        fileTuple
        sampleName


    main:
        infq = []
        infq = fileTuple.collect { ' --in-fq ' + it.fastq1 + ' ' + it.fastq2 + ' "' + it.rg + '" ' }

        pb_germline(infq.join(' '), sampleName)
    emit:
        pb_germline.out
}

workflow PB_DEEPVARIANT {
    take:         
        bamFile
        sampleName

    main:
        pb_deepvariant(bamFile, sampleName)
    emit:
        pb_deepvariant.out
}

workflow PB_HAPLOTYPECALLER {
    take:         
        bamFile
        sampleName

    main:
        pb_haplotypecaller(bamFile, sampleName)
    emit:
        pb_haplotypecaller.out
}

