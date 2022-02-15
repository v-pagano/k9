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
        val sampleName

    output:
        file "${sampleName}*"

    cpus params.samtoolsCpus
    container params.samtoolsContainer

    script:
    """
        samtools merge --threads ${params.samtoolsCpus} -c -f -l 6 '${sampleName}.bam' ${f.join(' ')}
        samtools index '${sampleName}.bam'
    """
}
