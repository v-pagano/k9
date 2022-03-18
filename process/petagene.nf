process bam2pgbam {
    
    input:
        tuple val(sample), path(bam)

    output:
        path "*.pgbam"
 
    cpus params.bam2pgbamCpus

    when:
        params.petagene

    script:
    """
        module load ${params.petageneModule}
        petasuite -c -m lossless -v off -s ${params.species} ${bam}
    """
}
