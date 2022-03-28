process sam_sort {
    input:
        path f
    output:
        file 'temp.bam'

    cpus params.samtoolsCpus
    container params.samtoolsContainer

    script:
        """
            samtools fixmate --threads ${params.samtoolsCpus} -m '${f}' - | \
            samtools sort -l 2 -m 3G --threads 24 --output-fmt BAM -o 'temp.bam'
        """

}

process sam_bam2fq {
    input:
        tuple val(sample), path(bam)

    output:
        tuple val(sample), path("${sample}_1.fq"), path("${sample}_2.fq")

    cpus params.samtoolsCpus
    container params.samtoolsContainer

    script:
    """
        samtools bam2fq --threads ${params.samtoolsCpus} -1 ${sample}_1.fq -2 ${sample}_2.fq '${bam}'
    """

}

process sam_merge {
    input:
        val f
        tuple val(sample), path(fq)
        val suffix

    output:
        tuple val(sample), path("${sample}_${suffix}.bam"), emit: bam
        path "${sample}*", emit: publishFiles

    cpus params.samtoolsCpus
    container params.samtoolsContainer

    script:
    """
        samtools merge --threads ${params.samtoolsCpus} -c -f -l 6 '${sample}_${suffix}.bam' ${f.join(' ')}
        samtools index '${sample}_${suffix}.bam'
    """
}
