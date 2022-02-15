import org.yaml.snakeyaml.Yaml

include { PB_HAPLOTYPECALLER; PB_DEEPVARIANT } from './parabricks'
include { haplotypecaller; vcf_merge; make_examples; call_variants; post_process_calls; filter_variants } from '../process/variantcallers'
include { ANNOTATION_NORMAL; ANNOTATION_TUMOR } from './annotation'

workflow VARIANTCALLERS_TUMOR {
    take:
        bamFile
        sampleName

    main:
        VARIANTCALLERS_NORMAL(bamFile, sampleName)
        ANNOTATION_TUMOR(VARIANTCALLERS_NORMAL.out, sampleName)

    emit:
        VARIANTCALLERS_NORMAL.out
}

workflow VARIANTCALLERS_NORMAL {
    take:
        bamFile
        sampleName

    main:
        vcfs = bamFile.flatten().filter{ it ==~ /.*vcf/ || it ==~ /.*vcf.gz/ }
        bams = bamFile.flatten().filter{ it ==~ /.*\.bam/ }
        // Handle accelerated parabricks germline
        if (!(params.accelerated || params.mapping_type == 'parabricks')) {
            if (params.pb_haplotypecaller) {
                PB_HAPLOTYPECALLER(bams, sampleName)
                vcfs = vcfs.concat(PB_HAPLOTYPECALLER.out.flatten().filter{ it ==~ /.*vcf/ || it ==~ /.*vcf.gz/ })
            }
        }
        PB_DEEPVARIANT(bams, sampleName)
        DEEPVARIANT(bams, sampleName, bamFile)
        vcfs = vcfs.concat(PB_DEEPVARIANT.out.flatten().filter{ it ==~ /.*vcf/ || it ==~ /.*vcf.gz/ })
        vcfs = vcfs.concat(DEEPVARIANT.out.flatten().filter{ it ==~ /.*vcf/ || it ==~ /.*vcf.gz/ })
        //Need this here because of the integrated germline pipeline for parabricks
        if (params.haplotypecaller) {
            HAPLOTYPECALLER(bams, sampleName)
            vcfs = vcfs.concat(HAPLOTYPECALLER.out.flatten().filter{ it ==~ /.*vcf/ || it ==~ /.*vcf.gz/ })
        }

        ANNOTATION_NORMAL(vcfs, sampleName)

    emit:
        vcfs
}

workflow HAPLOTYPECALLER {
    take: 
        bamFile
        sampleName

    main:
        Yaml parser = new Yaml()
        intervals = parser.load((params.HaplotypecallerIntervalsYaml as File).text)

        txtInt = ''
        for (interval in intervals.calling_intervals) {
            arrInt = interval.collect { ' -L "' + it.contig + ':' + it.start + '-' + it.stop + '" '}
        }
        
        haplotypecaller(bamFile, sampleName, arrInt)
        vcf_merge(haplotypecaller.out.collect(), sampleName)

    emit:
        vcf_merge.out
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