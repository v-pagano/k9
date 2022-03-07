process fastqc {
    
    container params.fastqcContainer

    input:
        tuple val(sample), path(fq)

    output:
        path "*.html"
        path "*.zip"
 
    when:
        params.fastqc

    cpus params.fastqcCPUs

    script:
    """
        fastqc --threads ${params.fastqcCPUs} ${fq[0]} ${fq[1]}
    """
}