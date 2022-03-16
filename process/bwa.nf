
process bwa_mem {

    input:
        tuple val(sample), path(fq)

    output:
        file "${sample}.bam"

    cpus params.bwamemCpus
    container params.bwaContainer

    when:
        params.bwa

    script:
        """
            bwa mem -v 3 -Y -K 100000000 -t ${params.bwamemCpus} \
            '${params.pb_reference}' \
            '${fq[0]}' ${fq[1]} \
            > '${sample}.bam'
        """

}

process bwa_mem2 {
    
    input:
        tuple val(sample), path(fq)

    output:
        file "${sample}.bam"

    cpus params.bwamemCpus
    container params.bwa2Container

    when:
        params.bwa2

    script:
        """
            bwa-mem2 mem -v 3 -Y -K 100000000 -t ${params.bwamemCpus} \
            '${params.bwa2reference}' \
            '${fq[0]}' \
            '${fq[1]}' \
            > '${sample}.bam'
        """

}
