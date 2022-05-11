process bam2pgbam {
    
    input:
        tuple val(sample), val(bam)

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

process petageneCompressBAM {

    input:
        tuple val(meta), path(bam)
        val species
        val encrypt
        val datasteward

    output:
        tuple val(meta), path("*.pgbam*"), emit: pgbam
        path("*.pgbam*"), emit: publishFiles

    cpus params.bam2pgbamCpus
    
    script:
    """
        source /etc/profile.d/modules.sh
        module load petagene/protect_1.3.11
        petasuite --compress ${encrypt ? '--encrypt' : ''} -m lossless -v off --datasteward ${datasteward} -s ${species} ${bam}
    """

}

process petageneCompressFastq {

    input:
        tuple val(meta), path(fq)
        val species
        val encrypt
        val datasteward

    output:
        tuple val(meta), path("*.fasterq*")

    cpus params.bam2pgbamCpus
    
    script:
    """
        source /etc/profile.d/modules.sh
        module load petagene/protect_1.3.11
        petasuite --compress ${encrypt ? '--encrypt' : ''} -m lossless -v off --datasteward ${datasteward} -s ${species} ${fq[0]}
        petasuite --compress ${encrypt ? '--encrypt' : ''} -m lossless -v off --datasteward ${datasteward} -s ${species} ${fq[1]}
    """

}

process petageneExpandFasterq {

    input:
        tuple val(meta), path(fq)
        val species
        val datasteward

    output:
        tuple val(meta), path("*.fastq.gz")
    
    cpus params.bam2pgbamCpus

    script:
    """
        source /etc/profile.d/modules.sh
        module load petagene/protect_1.3.11
        petasuite --decompress -m lossless -v off --datasteward ${datasteward} -s ${species} ${fq[0]}
        petasuite --decompress -m lossless -v off --datasteward ${datasteward} -s ${species} ${fq[1]}
    """

}

process petageneEncryptVCF {

    input:
        tuple val(meta), path(vcf)
        val species
        val encrypt
        val datasteward

    output:
        tuple val(meta), path("*_p"), emit: pgvcf
        path("*_p"), emit: publishFiles

    cpus params.bam2pgbamCpus
    
    script:
    """
        source /etc/profile.d/modules.sh
        module load petagene/protect_1.3.11
        petasuite --encrypt --datasteward ${datasteward} ${vcf}
    """

}
