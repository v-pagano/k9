process fastqc {
    
    container params.fastqcContainer

    input:
        tuple val(sample), path(fq)

    output:
        path "*.html"
        path "*.zip"
 
    cpus params.fastqcCPUs

    when:
        params.fastqc


    script:
    """
        fastqc --threads ${params.fastqcCPUs} ${fq[0]} ${fq[1]}
    """
}