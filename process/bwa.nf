
process bwa_mem {

    input:
        tuple val(sample), path(fq)

    output:
        file "${sample}_bwa.bam"

    cpus params.bwamemCpus
    container params.bwaContainer

    when:
        params.bwa

    script:
        """
            ${params.petagene ? 'LD_PRELOAD=' + params.petalinkModule : ''} \
            bwa mem -v 3 -Y -K 100000000 -t ${params.bwamemCpus} \
            '${params.pb_reference}' \
            -R '@RG\\tID:${sample}\\tLB:lib1\\tPL:bar\\tSM:sample\\tPU:${sample}' \
            '${fq[0]}' \
            '${fq[1]}' \
            > '${sample}_bwa.bam'
        """

}

process bwa_mem2 {
    
    input:
        tuple val(sample), path(fq)

    output:
        file "${sample}_bwa2.bam"

    cpus params.bwamemCpus
    container params.bwa2Container

    when:
        params.bwa2

    script:
        """
            bwa-mem2 mem -v 3 -Y -K 100000000 -t ${params.bwamemCpus} \
            '${params.bwa2reference}' \
            -R '@RG\\tID:${sample}\\tLB:lib1\\tPL:bar\\tSM:sample\\tPU:${sample}' \
            '${fq[0]}' \
            '${fq[1]}' \
            > '${sample}_bwa2.bam'
        """

}
