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

process trimmomatic {

    script:
    """
        java -jar $TRIMMOMATIC_PATH/trimmomatic-0.33.jar PE -threads 1 -phred33 $reads $prefix"_filtered_R1.fastq" $prefix"_unpaired_R1.fastq" $prefix"_filtered_R2.fastq" $prefix"_unpaired_R2.fastq" ILLUMINACLIP:${params.trimmomatic_adapters_file}:${params.trimmomatic_adapters_parameters} SLIDINGWINDOW:${params.trimmomatic_window_length}:${params.trimmomatic_window_value} MINLEN:${params.trimmomatic_mininum_length} 2> ${name}.log
    """
}