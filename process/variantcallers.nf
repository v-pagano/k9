process haplotypecaller {
    input:
        tuple val(sample), path(bam)
        each interval

    output:
        file "${sample}.g.vcf.gz"

    container params.gatkContainer

    when:
        params.haplotypecaller

    cpus params.haplotypecallerCpus

    script:
        """
            gatk HaplotypeCaller \
            --input "${sample}.bam" \
            -O "${sample}.g.vcf.gz" \
            --reference ${params.reference} \
            --java-options "-Xmx14G" \
            -ERC GVCF \
            ${interval}
        """
}

process vcf_merge {
    input:
        val vcfList
        val sample

    output:
        file "${sample}*"

    container params.gatkContainer

    script:
    vcfList.view()
    """
        gatk MergeVcfs \
        --java-options "-Xmx7G" \
        --INPUT ${vcfList.join(' --INPUT ')} \
        --OUTPUT "${sample}_haplotypecaller.g.vcf.gz"
    """
}

process make_examples {
    input:
        val bamFile
        val sampleName
        each instance
        file extraFiles

    output:
        file "*.gz"

    cpus params.make_examplesCpus
    container params.deepvariantContainer


    when:
        params.deepvariant

    script:
    """
        make_examples \
        --logging_level WARN \
        --mode calling \
        --task ${instance} \
        --ref "${params.reference}" \
        --reads "${bamFile}" \
        --examples "${sampleName}.ex.tfrecord@${params.make_examplesInstances}.gz" \
        --gvcf "${sampleName}.gvcf.tfrecord@${params.make_examplesInstances}.gz"
    """
   
}

process call_variants {

    input:
        path exampleFiles
        val sampleName        

    output:
        file "*.gz"

    queue params.gpuPartition
    cpus params.deepvariantCpus
    clusterOptions '--gres gpu:1 -N 1 --tasks-per-node 1'
    container params.deepvariantContainer

    when:
        params.deepvariant

    script:
    """
        call_variants \
            --checkpoint  "${params.deepvariantModel}" \
            --examples "${sampleName}.ex.tfrecord@${params.make_examplesInstances}.gz" \
            --outfile "${sampleName}.cvo.tfrecord.gz"
    """
}

process post_process_calls {

    input:
        path exampleFiles
        path variantFiles
        val sampleName        

    output:
        file "*.gz"

    cpus params.deepvariantPostProcessCpus
    container params.deepvariantContainer

    when:
        params.deepvariant

    script:
    """
        postprocess_variants \
            --ref "${params.reference}" \
            --infile "${sampleName}.cvo.tfrecord.gz" \
            --outfile "${sampleName}_deepvariant.vcf.gz" \
            --nonvariant_site_tfrecord_path "${sampleName}.gvcf.tfrecord@${params.make_examplesInstances}.gz" \
            --gvcf_outfile "${sampleName}_deepvariant.g.vcf.gz"

        bcftools index --tbi --force "${sampleName}_deepvariant.vcf.gz"
    """
}

process filter_variants {

    input:
        path variantFiles
        val sampleName        

    output:
        tuple path("${sampleName}_deepvariant.vcf.gz"), path("${sampleName}_deepvariant_pass.vcf.gz"), path("${sampleName}_deepvariant.g.vcf.gz")

    cpus params.deepvariantPostProcessCpus
    container params.deepvariantContainer

    when:
        params.deepvariant

    script:
    """
        bcftools filter \
            --output-type z \
            --include 'FILTER == "PASS"' \
            "${sampleName}_deepvariant.vcf.gz" \
            > "${sampleName}_deepvariant_pass.vcf.gz"

        bcftools index --tbi --force "${sampleName}_deepvariant_pass.vcf.gz"
    """
}