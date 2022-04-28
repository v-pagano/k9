import org.yaml.snakeyaml.Yaml

include { pb_haplotypecaller; pb_deepvariant; pb_cnvkit; pb_manta_normalOnly; pb_strelka_normalOnly  } from '../process/parabricks'
include { haplotypecaller; vcf_merge; make_examples; call_variants; post_process_calls; filter_variants } from '../process/variantcallers'

workflow VARIANTCALLERS {
    take:
        bam

    main:

        vcfFiles = Channel.empty()
        publishFiles = Channel.empty()

        pb_haplotypecaller(bam, params.pb_reference)
        publishFiles = publishFiles.mix(pb_haplotypecaller.out[1].flatten())
        vcfFiles = vcfFiles.mix(pb_haplotypecaller.out[0])
        
        pb_deepvariant(bam, params.pb_reference)
        publishFiles = publishFiles.mix(pb_deepvariant.out[1].flatten())
        vcfFiles = vcfFiles.mix(pb_deepvariant.out[0])

        pb_cnvkit(bam, params.pb_reference)
        publishFiles = publishFiles.mix(pb_cnvkit.out.flatten())

        pb_manta_normalOnly(bam, params.pb_reference)
        publishFiles = publishFiles.mix(pb_manta_normalOnly.out[1].flatten())

        pb_strelka_normalOnly(bam, params.pb_reference)
        publishFiles = publishFiles.mix(pb_strelka_normalOnly.out[1].flatten())

        if (params.haplotypecaller) {
            HAPLOTYPECALLER(bam)
            publishFiles = publishFiles.mix(HAPLOTYPECALLER.out.flatten())
            vcfFiles = vcfFiles.mix(HAPLOTYPECALLER.out)
        }

    emit:
        vcf = vcfFiles
        publishFiles = publishFiles
}

workflow HAPLOTYPECALLER {
    take: 
        bam

    main:

        Yaml parser = new Yaml()
        intervals = parser.load((params.HaplotypecallerIntervalsYaml as File).text)

        txtInt = ''
        arrInt = []
        // for (interval in intervals.calling_intervals) {
        //     for (i in interval) {
        //         txtInt = ' -L "' + i.contig + ':' + i.start + '-' + i.stop + '" '
        //         arrInt.add(txtInt)
        //     }
        // }

        for (interval in intervals.calling_intervals) {
            txtInt = interval.collect { ' -L "' + it.contig + ':' + it.start + '-' + it.stop + '" '}.join(' ')
            arrInt.add(txtInt)
        }

        haplotypecaller(bam, arrInt)
        sample = bam.map { it[0] }
        vcfFiles = haplotypecaller.out.map { it[1] }
        tempVCF = vcfFiles.collect()
        output = vcf_merge(sample, tempVCF)

    emit:
        output
}

workflow DEEPVARIANT {
    take: 
        bamFile
        sampleName
        extraFiles

    main:
        examples = make_examples(bamFile, sampleName, Channel.of(0..params.make_examplesInstances - 1), extraFiles)
        variants = call_variants(examples.collect(), sampleName)
        processedVariants = post_process_calls(examples.collect(), variants, sampleName)
        filter_variants(processedVariants, sampleName)

    emit:
        filter_variants.out
}