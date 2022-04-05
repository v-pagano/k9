nextflow.enable.dsl=2

include { MAPPING } from './workflows/mapping'
include { VARIANTCALLERS } from './workflows/variantcallers'
include { FQ_PREP } from './workflows/fq_prep'
include { BAMQC } from './workflows/bamqc'
include { ANNOTATION } from './workflows/annotation'
include { PB_HAPLOTYPECALLER; PB_DEEPVARIANT; PB_SOMATIC } from './workflows/parabricks'
include { BAM2PGBAM } from './workflows/petagene'
include { UBAM2FQ; FQ2UBAM } from './workflows/ubam'
include { publishResults; dot2svg } from './process/helpers'
include { sam_bam2fq } from './process/sam'

workflow {

    if (params.usegpu03) {
        params.gpuClusterOptions = '--nodelist=dback-gpu03'
        params.gpuPartition = 'gpu-dev'
    }
    
    publishFiles = Channel.empty()
    bamFiles = Channel.empty()
    vcfFiles = Channel.empty()

    if (params.inputType == 'vcf') {

        vcfFiles = Channel.fromPath(params.input)
        vcfs = vcfFiles.map { [ it.baseName, it ] }
    
        ANNOTATION(vcfs)
        publishFiles = publishFiles.mix(ANNOTATION.out)

        publishResults(publishFiles.flatten())

    }

    if (params.inputType == 'bam') {

        bamFiles = Channel.fromPath(params.input)
        bams = bamFiles.map { [ it.baseName, it ] }

        BAMQC(bams)
        publishFiles = publishFiles.mix(BAMQC.out)

        VARIANTCALLERS(bams)
        publishFiles = publishFiles.mix(VARIANTCALLERS.out[1])

        ANNOTATION(VARIANTCALLERS.out[0])
        publishFiles = publishFiles.mix(ANNOTATION.out)

        BAM2PGBAM(bams)
        publishFiles = publishFiles.mix(BAM2PGBAM.out)

        publishResults(publishFiles.flatten())

    }

    if (params.inputType == 'ubam') {

        bamFiles = Channel.fromPath(params.input)
        bams = bamFiles.map { [ it.baseName, it ] }

        UBAM2FQ(bams)
        fastq = UBAM2FQ.out.map { [ it[0] , [ it[1], it[2] ] ] }

    } 
    
    if (params.inputType == 'fastq') {

        fastq = Channel.fromFilePairs(params.input)

    }

    if (params.inputType == 'fastq' || params.inputType == 'ubam') {

        FQ_PREP(fastq)
        publishFiles = publishFiles.mix(FQ_PREP.out)

        FQ2UBAM(fastq)
        publishFiles = publishFiles.mix(FQ2UBAM.out)

        MAPPING(fastq)
        publishFiles = publishFiles.mix(MAPPING.out.publishFiles)

        BAMQC(MAPPING.out.bam)
        publishFiles = publishFiles.mix(BAMQC.out)

        VARIANTCALLERS(MAPPING.out.bam)
        publishFiles = publishFiles.mix(VARIANTCALLERS.out[1])

        ANNOTATION(VARIANTCALLERS.out[0])
        publishFiles = publishFiles.mix(ANNOTATION.out)

        BAM2PGBAM(MAPPING.out.bam)
        publishFiles = publishFiles.mix(BAM2PGBAM.out)

        publishResults(publishFiles.flatten())

    }

}