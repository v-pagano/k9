include { PB_FQ2BAM; PB_GERMLINE } from './parabricks'
include { BWA; BWA2 } from './bwa'

workflow MAPPING_TUMOR {
    take:
        sample

    main:
        MAPPING_NORMAL(sample)

    emit:
        MAPPING_NORMAL.out
}

workflow MAPPING_NORMAL {
    take:
        sample

    main:
        if (params.accelerated || params.mapping_type == 'parabricks') {
            if (params.pb_haplotypecaller) {
                bamOutput = PB_GERMLINE(sample, sample[0].sample)
            } else {
                bamOutput = PB_FQ2BAM(sample, sample[0].sample)
            }
        } else {
            switch(params.mapping_type) {
                case 'bwa':
                    bamOutput = BWA(sample, sample[0].sample)
                    break
                case 'bwa2':
                    bamOutput = BWA2(sample, sample[0].sample)
                    break
            }
        }

    emit:
        bamOutput
}
