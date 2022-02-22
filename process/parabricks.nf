process pb_fq2bam {
    input:
        val fq
        val sampleName

    output:
        file "${sampleName}*"

    queue params.gpuPartition
    clusterOptions '--exclusive'

    script:
    """
        source /etc/profile.d/modules.sh
        module load parabricks/${params.pb_ver} 
        pbrun fq2bam --bwa-options '-K 100000000 -Y' --ref ${params.pb_reference} \
        ${fq} \
        --out-bam '${sampleName}.bam' \
        ${params.baserecalibration ? '--knownSites ' + params.knownSites + ' --out-recal-file ' + sampleName + '_recal.txt ' : ''} --tmp-dir /scratch/vpagano/tmp
    """

}

process pb_rna_fq2bam {
    input:
        val f
    output:
        tuple file('merged.bam'), file('results/*')

    queue params.gpuPartition
    clusterOptions '--exclusive'

    script:
    """
        source /etc/profile.d/modules.sh
        module load parabricks/${pb_ver}
        mkdir results
        pbrun rna_fq2bam --ref ${params.pb_reference} \
        ${f.join(' ')} \
        --read-files-command zcat \
        --genome-lib-dir /home/tgenref/homo_sapiens/grch38_hg38/hg38tgen/gene_model/ensembl_v98/tool_resources/star_2.7.3a/83bpReads \
        --sjdb-overhang 82 \
        --out-bam 'merged.bam' \
        --output-dir results \
        --tmp-dir /scratch/vpagano/tmp
    """

}

process pb_deepvariant {
    input:
        file bam
        val sampleName
    output:
        file "${sampleName}*"

    queue params.gpuPartition
    clusterOptions '--exclusive'

    when:
        params.pb_deepvariant

    script:
    """
        source /etc/profile.d/modules.sh
        module load parabricks/${params.pb_ver} 
        pbrun deepvariant --ref ${params.pb_reference} \
        --in-bam '${sampleName}.bam' \
        --out-variants '${sampleName}_pb_deepvariant.vcf' \
        --tmp-dir /scratch/vpagano/tmp
    """

}

process pb_haplotypecaller {
    input:
        file bam
        val sampleName
    output:
        file "${sampleName}*"

    queue params.gpuPartition
    clusterOptions '--exclusive'

    script:
    """
        source /etc/profile.d/modules.sh
        module load parabricks/${params.pb_ver} 
        pbrun haplotypecaller --ref ${params.pb_reference} \
        --in-bam '${sampleName}.bam' \
        --out-variants '${sampleName}_pb_haplotypecaller.vcf' \
        --tmp-dir /scratch/vpagano/tmp
    """

}

process pb_germline {
    input:
        each fq
    output:
        file "${fq[0]}*"

    queue params.gpuPartition
    clusterOptions '--exclusive'
    publishDir '/scratch/vpagano/results/canine', mode: 'copy', overwrite: 'true'

    script:
    """
        source /etc/profile.d/modules.sh
        module load parabricks/${params.pb_ver} 
        pbrun germline --bwa-options '-K 100000000 -Y' --ref ${params.pb_reference} \
        --in-fq ${fq[1][0]} ${fq[1][1]} \
        --out-bam '${fq[0]}.bam' \
        --out-variants '${fq[0]}_pb_haplotypecaller.vcf' \
        --tmp-dir /scratch/vpagano/tmp
    """

}

process pb_somatic {
    input:
        val fqTumor
        val fqNormal
        val sampleNameTumor
        val sampleNameNormal
    output:
        tuple file("${sampleNameTumor}*"), file("${sampleNameNormal}*")
        
    queue params.gpuPartition
    clusterOptions '--exclusive'

    script:
    """
        source /etc/profile.d/modules.sh
        module load parabricks/${params.pb_ver} 
        pbrun somatic --ref ${params.pb_reference} \
        --in-tumor-fq ${fqTumor} \
        --out-vcf '${sampleNameTumor}_${sampleNameNormal}.vcf' \
        --out-tumor-bam '${sampleNameTumor}.bam' \
        --in-normal-fq ${fqNormal} \
        --out-normal-bam '${sampleNameNormal}.bam'  \
        --tmp-dir /scratch/vpagano/tmp
    """

}

process mutect {

    input:
        file tumor_bam
        file normal_bam

    output:

    queue params.gpuPartition
    clusterOptions '--exclusive'

    script:
    """
        pbrun mutectcaller \
        --ref ${params.reference} \
        --in-tumor-bam ${tumor_bam} \
        --tumor-name {{ pair.tumor.rgsm }} \
        --in-normal-bam ${normal_bam} \
        --normal-name {{ pair.normal.rgsm }} \
        --tmp-dir /scratch/vpagano/tmp \
        --out-vcf mutect.vcf
    """
}

