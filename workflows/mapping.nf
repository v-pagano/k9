include { pb_fq2bam; pb_germline; pb_somaticTumorOnly } from '../process/parabricks'
include { BWA; BWA2 } from './bwa'
include { PB_FQ2BAM; PB_GERMLINE; PB_DEEPVARIANT_GERMLINE } from './parabricks'

workflow MAPPING {
    take:
        fastq

    main:

        publishFiles = Channel.empty()
        bamFiles = Channel.empty()
        vcfFiles = Channel.empty()

        if (params.pb_somatic) {
            if (params.somaticType = 'tumorOnly') {
                pb_somaticTumorOnly(fastq, params.pb_reference)
                publishFiles = publishFiles.mix(pb_somaticTumorOnly.out.publishFiles.flatten())
                bamFiles = bamFiles.mix(pb_somaticTumorOnly.out.bam)
                vcfFiles = vcfFiles.mix(pb_somaticTumorOnly.out.vcf)
            } else {

            }
        } else {

            PB_FQ2BAM(fastq, params.pb_reference) 
            publishFiles = publishFiles.mix(PB_FQ2BAM.out.publishFiles.flatten())
            bamFiles = bamFiles.mix(PB_FQ2BAM.out.bamFiles)

            PB_GERMLINE(fastq, params.pb_reference)
            publishFiles = publishFiles.mix(PB_GERMLINE.out.publishFiles.flatten())
            bamFiles = bamFiles.mix(PB_GERMLINE.out.bam)
            vcfFiles = vcfFiles.mix(PB_GERMLINE.out.vcf)

            PB_DEEPVARIANT_GERMLINE(fastq, params.pb_reference)
            publishFiles = publishFiles.mix(PB_DEEPVARIANT_GERMLINE.out.publishFiles.flatten())
            bamFiles = bamFiles.mix(PB_DEEPVARIANT_GERMLINE.out.bam)
            vcfFiles = vcfFiles.mix(PB_DEEPVARIANT_GERMLINE.out.vcf)

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
        }

    emit:
        bam = bamFiles
        publishFiles = publishFiles
        vcf = vcfFiles
}
