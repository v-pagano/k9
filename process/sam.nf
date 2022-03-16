process sam_sort {
    input:
        each f
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

process sam_merge {
    input:
        val f
        tuple val(sample), path(fq)

    output:
        file "${sample}*"

    cpus params.samtoolsCpus
    container params.samtoolsContainer

    script:
    """
        samtools merge --threads ${params.samtoolsCpus} -c -f -l 6 '${sample}.bam' ${f.join(' ')}
        samtools index '${sample}.bam'
    """
}
