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

        if (params.bwa) {
            BWA(fastq)
            publishFiles = publishFiles.mix(BWA.out.publishFiles.flatten())
            bamFiles = bamFiles.mix(BWA.out.bam)
        }
                
        if (params.bwa2) {
            BWA2(fastq)
            publishFiles = publishFiles.mix(BWA2.out[1].flatten())
            bamFiles = bamFiles.mix(BWA2.out[0])
        }
        
    emit:
        bam = bamFiles
        publishFiles = publishFiles
        vcf = vcfFiles
}
