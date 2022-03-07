process vep {

    input:
        tuple val(sample), path(vcf)
        path reference

    output:
        file "${sample}*"

    cpus params.vepCpus
    container params.vepContainer

    when:
        params.vep

    script:
    //Split up these processes when containers are ready
    """
        vep --fork 4 \
        --input_file "${vcf}" \
        --format vcf \
        --output_file "${sample}_vep.vcf" \
        --vcf \
        --vcf_info_field CSQ \
        --species homo_sapiens \
        --force_overwrite \
        --no_stats \
        --cache \
        --dir_cache "${params.vepData}" \
        --cache_version 98 \
        --offline \
        --fasta "${reference}" \
        --buffer_size 10000 \
        --terms SO \
        --hgvs \
        --hgvsg \
        --symbol \
        --sift b \
        --polyphen b \
        --humdiv \
        --uniprot \
        --domains \
        --canonical \
        --flag_pick_allele_gene \
        --pick_order canonical,appris,tsl,biotype,rank,ccds,length

        bcftools view \
        --threads 4 \
        --output-type z \
        --output-file "${sample}_vep.vcf.gz" \
        "${sample}_vep.vcf"

        rm "${sample}_vep.vcf"

        bcftools index --threads 4 --force --tbi "${sample}_vep.vcf.gz"
    """
}

process snpeff {

    input:
        tuple val(sample), path(vcf)

    output:
        file "${sample}*"

    cpus params.snpeffCpus
    container params.snpeffContainer

    when:
        params.snpeff

    script:
    //Split up these processes when containers are ready
    """
        snpEff ann \
            -t \
            -c "${params.snpeffConfig}" \
            -dataDir "${params.snpeffData}" \
            -canon \
            -hgvs \
            -lof \
            "${params.snpeffDb}" \
            "${vcf}" \
            | bcftools view --output-type z --output-file "${sample}_snpeff.vcf.gz"

        bcftools index --threads 4 --tbi --force "${sample}_snpeff.vcf.gz"

    """
}