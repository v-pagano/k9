include { pb_fq2bam; pb_germline } from '../process/parabricks'
include { BWA; BWA2 } from './bwa'

workflow MAPPING {
    take:
        fastq

    main:

        publishFiles = Channel.empty()
        bamFiles = Channel.empty()
        vcfFiles = Channel.empty()

        pb_fq2bam(fastq, params.pb_reference) 
        publishFiles = publishFiles.mix(pb_fq2bam.out[1].flatten())
        bamFiles = bamFiles.mix(pb_fq2bam.out[0])

        pb_germline(fastq, params.pb_reference)
        publishFiles = publishFiles.mix(pb_germline.out[1].flatten())
        bamFiles = bamFiles.mix(pb_germline.out[0])
        vcfFiles = vcfFiles.mix(pb_germline.out[2])
                
        // pb_haplotypecaller(pb_fq2bam.out[0], params.pb_reference)
        // publishFiles = publishFiles.mix(pb_haplotypecaller.out[1].flatten())
        // vcfFiles = vcfFiles.mix(pb_haplotypecaller.out[0])
        
        // pb_deepvariant(pb_fq2bam.out[0], params.pb_reference)
        // publishFiles = publishFiles.mix(pb_deepvariant.out[1].flatten())
        // vcfFiles = vcfFiles.mix(pb_deepvariant.out[0])

    emit:
        bamFiles
        publishFiles
        vcfFiles
}
