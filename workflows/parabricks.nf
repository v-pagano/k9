include { pb_fq2bam; pb_germline; pb_deepvariant_germline; pb_haplotypecaller; pb_deepvariant; pb_somatic } from '../process/parabricks.nf'
include { petageneCompressBAM; petageneEncryptVCF } from '../process/petagene.nf'
include { gzip } from '../process/helpers.nf'

workflow PB_FQ2BAM {
    take: 
        fq
        reference

    main:
        pb_fq2bam(fq, reference)
        bamfiles = pb_fq2bam.out.bam
        if (params.petagene) {
            petageneCompressBAM(pb_fq2bam.out.bam, params.species, params.encrypt, params.datasteward)
            publishFiles = petageneCompressBAM.out.publishFiles
        } else {
            publishFiles = pb_fq2bam.out.publishFiles
        }

    emit:
        bamFiles = bamfiles
        publishFiles = publishFiles
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
        fq
        reference

    main:
        pb_germline(fq, reference)
        bamfiles = pb_germline.out.bam
        if (params.petagene) {
            petageneCompressBAM(pb_germline.out.bam, params.species, params.encrypt, params.datasteward)
            petageneEncryptVCF(pb_germline.out.vcf, params.species, params.encrypt, params.datasteward)
            publishFiles = petageneCompressBAM.out.publishFiles.mix(petageneEncryptVCF.out.publishFiles)
            vcfout = pb_germline.out.vcf
        } else {
            gzip(pb_germline.out.vcf)
            publishFiles = pb_germline.out.publishFiles.mix(gzip.out)
            vcfout = gzip.out
        }
        
    emit:
        bam = bamfiles
        publishFiles = publishFiles
        vcf = vcfout
}

workflow PB_DEEPVARIANT_GERMLINE {
    take: 
        fq
        reference

    main:
        pb_deepvariant_germline(fq, reference)
        bamfiles = pb_deepvariant_germline.out.bam
        if (params.petagene) {
            petageneCompressBAM(pb_deepvariant_germline.out.bam, params.species, params.encrypt, params.datasteward)
            petageneEncryptVCF(pb_deepvariant_germline.out.vcf, params.species, params.encrypt, params.datasteward)
            publishFiles = petageneCompressBAM.out.publishFiles.mix(petageneEncryptVCF.out.publishFiles)
            vcfout = pb_deepvariant_germline.out.vcf
        } else {
            gzip(pb_deepvariant_germline.out.vcf)
            publishFiles = pb_deepvariant_germline.out.publishFiles.mix(gzip.out)
            vcfout = gzip.out
        }
        
    emit:
        bam = bamfiles
        publishFiles = publishFiles
        vcf = vcfout
}


workflow PB_DEEPVARIANT {
    take:         
        bamFile
        sampleName

    main:
        pb_deepvariant(bamFile)

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

