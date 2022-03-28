process picard_fq2ubam {

    input:
        tuple val(sample), path(fq)

    output:
        tuple val(sample), path("${sample}.ubam")

    container params.picardContainer

    cpus params.samtoolsCpus

    when:
        params.ubam

    script:
    """
        java -Xmx8G -jar \
        /usr/local/share/picard-2.26.11-0/picard.jar \
        FastqToSam \
        FASTQ=${fq[0]} \
        FASTQ2=${fq[1]} \
        OUTPUT=${sample}.ubam \
        READ_GROUP_NAME=H0164.2 \
        SAMPLE_NAME=NA12878 \
        LIBRARY_NAME=Solexa-272222 \
        PLATFORM_UNIT=H0164ALXX140820.2 \
        PLATFORM=illumina \
        SEQUENCING_CENTER=BI
    """

}
