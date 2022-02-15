
process bwa_mem {
    input:
        each f
    output:
        file "${f.sample}.bam"

    cpus params.bwamemCpus
    container params.bwaContainer

    script:
        """
            bwa mem -v 3 -Y -K 100000000 -t ${params.bwamemCpus} -R '${f.rg}' \
            '${params.reference}' \
            '${f.fastq1}' ${f.fastq2} \
            > '${f.sample}.bam'
        """

}

process bwa_mem2 {
    input:
        each f
    output:
        file "${f.sample}.bam"

    cpus params.bwamemCpus
    container params.bwa2Container

    script:
        """
            bwa-mem2 mem -v 3 -Y -K 100000000 -t ${params.bwamemCpus} -R '${f.rg}' \
            '${params.bwa2reference}' \
            '${f.fastq1}' \
            '${f.fastq2}' \
            > '${f.sample}.bam'
        """

}
