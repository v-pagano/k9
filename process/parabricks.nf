process pb_fq2bam {
    input:
        tuple val(sample), path(fq)
        path reference

    output:
        tuple val("${sample}"), path("${sample}_pb.bam"), emit: bam
        path("${sample}*"), emit: publishFiles

    queue (params.usegpu03 ? 'gpu-dev' : params.gpuPartition)
    clusterOptions (params.usegpu03 ? '--nodelist=dback-gpu03 --exclusive' : "--exclusive ${params.gpuClusterOptions}")

    when:
        params.parabricks && !params.pb_haplotypecaller && !params.pb_deepvariant 

    script:
    """
        source /etc/profile.d/modules.sh
        module load parabricks/${params.pb_ver}
        ${params.petagene ? 'LD_PRELOAD=' + params.petalinkModule : ''} \
        ${params.petagene ? 'SINGULARITYENV_LD_PRELOAD=' + params.petalinkModule : ''} \
        pbrun fq2bam --bwa-options '-K 100000000 -Y' --ref ${reference} \
        --in-fq ${fq[0]} ${fq[1]} \
        --out-bam '${sample}_pb.bam' \
        ${params.petagene ? '--with-petagene-dir ' + params.petageneFolder : ''} \
        ${params.baserecalibration ? '--knownSites ' + params.knownSites + ' --out-recal-file ' + sampleName + '_recal.txt ' : ''} \
        --tmp-dir /scratch/vpagano/tmp
    """

}

process pb_rna_fq2bam {
    input:
        val f
    output:
        tuple file('merged.bam'), file('results/*')

    queue params.gpuPartition
    clusterOptions "--exclusive ${params.gpuClusterOptions}"

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
        tuple val(sample), path(bam)
        path reference

    output:
        tuple val("${sample}"), path("${sample}*deepvariant.vcf")
        file "${sample}*"

    queue params.gpuPartition
    clusterOptions "--exclusive ${params.gpuClusterOptions}"

    when:
        (params.pb_deepvariant && !params.parabricks && params.inputType == 'fastq') || (params.inputType == 'bam' && params.pb_deepvariant)

    script:
    """
        source /etc/profile.d/modules.sh
        module load parabricks/${params.pb_ver} 
        pbrun deepvariant --ref ${reference} \
        --in-bam '${bam}' \
        ${params.gvcf ? '--gvcf ' : ' '} \
        --out-variants '${sample}_pb_deepvariant.${params.gvcf ? 'g.' : ''}vcf' \
        --tmp-dir /scratch/vpagano/tmp
    """

}

process pb_haplotypecaller {
    input:
        tuple val(sample), path(bam)
        path reference

    output:
        tuple val("${sample}"), path("${sample}*.vcf")
        file "${sample}*"

    queue params.gpuPartition
    clusterOptions "--exclusive ${params.gpuClusterOptions}"

    //This boolean is complicated because we only want to run haplotypecaller if we are not aligning fastqs as well
    //Otherwise pb_germline is faster
    when:
        (params.pb_haplotypecaller && !params.parabricks && params.inputType == 'fastq') || (params.inputType == 'bam' && params.pb_haplotypecaller)

    script:
    """
        source /etc/profile.d/modules.sh
        module load parabricks/${params.pb_ver} 
        pbrun haplotypecaller --ref ${reference} \
        --in-bam '${bam}' \
        ${params.gvcf ? '--gvcf ' : ' '} \
        --out-variants '${sample}_pb_haplotypecaller.${params.gvcf ? 'g.' : ''}vcf' \
        --tmp-dir /scratch/vpagano/tmp
    """

}

process pb_germline {
    input:
        tuple val(sample), val(fq)
        path reference

    output:
        tuple val("${sample}"), path("${sample}_pb.bam"), emit: bam
        path("${sample}*"), emit: publishFiles
        tuple val("${sample}"), path("${sample}*.vcf"), emit: vcf

    queue (params.usegpu03 ? 'gpu-dev' : params.gpuPartition)
    clusterOptions (params.usegpu03 ? '--nodelist=dback-gpu03 --exclusive' : "--exclusive ${params.gpuClusterOptions}")

    when:
        params.parabricks && params.pb_haplotypecaller

    script:
    """
        source /etc/profile.d/modules.sh
        module load parabricks/${params.pb_ver} 
        ${params.petagene ? 'LD_PRELOAD=' + params.petalinkModule : ''} \
        ${params.petagene ? 'SINGULARITYENV_LD_PRELOAD=' + params.petalinkModule : ''} \
        pbrun germline --bwa-options '-K 100000000 -Y' --ref ${reference} \
        --in-fq ${fq[0]} ${fq[1]} \
        --out-bam '${sample}_pb.bam' \
        ${params.gvcf ? '--gvcf ' : ' '} \
        ${params.petagene ? '--with-petagene-dir ' + params.petageneFolder : ''} \
        --out-variants '${sample}_pb_haplotypecaller.${params.gvcf ? 'g.' : ''}vcf' \
        --tmp-dir /scratch/vpagano/tmp
    """

}

process pb_deepvariant_germline {
    input:
        tuple val(sample), val(fq)
        path reference

    output:
        tuple val("${sample}"), path("${sample}_pb.bam"), emit: bam
        path("${sample}*"), emit: publishFiles
        tuple val("${sample}"), path("${sample}*.vcf"), emit: vcf

    queue (params.usegpu03 ? 'gpu-dev' : params.gpuPartition)
    clusterOptions (params.usegpu03 ? '--nodelist=dback-gpu03 --exclusive' : "--exclusive ${params.gpuClusterOptions}")

    when:
        params.parabricks && params.pb_deepvariant

    script:
    """
        source /etc/profile.d/modules.sh
        module load parabricks/${params.pb_ver} 
        pbrun germline --bwa-options '-K 100000000 -Y' --ref ${reference} \
        --in-fq ${fq[0]} ${fq[1]} \
        --out-bam '${sample}_pb.bam' \
        ${params.gvcf ? '--gvcf ' : ' '} \
        --out-variants '${sample}_pb_deepvariant.${params.gvcf ? 'g.' : ''}vcf' \
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
    clusterOptions "--exclusive ${params.gpuClusterOptions}"

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
    clusterOptions "--exclusive ${params.gpuClusterOptions}"

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

process pb_collectmetrics {

    input:
        tuple val(sample), path(bam)
        path reference

    output:
        file "${sample}-qc/"

    queue params.gpuPartition
    clusterOptions "--exclusive ${params.gpuClusterOptions}"

    when:
        params.pb_collectmetrics

    script:
    """
        source /etc/profile.d/modules.sh
        module load parabricks/${params.pb_ver} 
        pbrun collectmultiplemetrics --ref ${reference} \
        --bam ${bam} \
        --out-qc-metrics-dir '${sample}-qc' \
        --gen-all-metrics \
        --tmp-dir /scratch/vpagano/tmp

    """
}

process pb_somaticTumorOnly {
    input:
        tuple val(sample), path(fq)
        path reference

    output:
        tuple val("${sample}"), path("${sample}_tumor_pb.bam"), emit: bam
        path("${sample}*"), emit: publishFiles
        tuple val("${sample}"), path("${sample}_pb_mutect.vcf"), emit: vcf
        
    queue (params.usegpu03 ? 'gpu-dev' : params.gpuPartition)
    clusterOptions (params.usegpu03 ? '--nodelist=dback-gpu03 --exclusive' : "--exclusive ${params.gpuClusterOptions}")

    script:
    """
        source /etc/profile.d/modules.sh
        module load parabricks/${params.pb_ver} 
        pbrun somatic --ref ${params.pb_reference} \
        --in-tumor-fq ${fq[0]} ${fq[1]} \
        --out-vcf '${sample}_pb_mutect.vcf' \
        --out-tumor-bam '${sample}_tumor_pb.bam' \
        --tmp-dir /scratch/vpagano/tmp
    """

}

process pb_cnvkit {

    input:
        tuple val(sample), path(bam)
        path reference

    output:
        file "${sample}-cnvkit/"

    cpus params.cnvkitCpus

    when:
        params.pb_cnvkit

    script:
    """
        source /etc/profile.d/modules.sh
        module load parabricks/${params.pb_ver} 
        pbrun cnvkit --ref ${reference} \
        --in-bam ${bam} \
        --output-dir '${sample}-cnvkit' \
        --tmp-dir /scratch/vpagano/tmp

    """
}

process pb_manta_normalOnly {
    input:
        tuple val(sample), path(bam)
        path reference

    output:
        tuple val("${sample}"), path("${sample}*")
        path "${sample}_pb_manta.manta_work"

    cpus params.mantaCpus

    when:
        params.pb_manta


    script:
    """
        source /etc/profile.d/modules.sh
        module load parabricks/${params.pb_ver} 
        pbrun manta --ref ${reference} \
        --in-normal-bam '${bam}' \
        --num-threads ${task.cpus} \
        --out-prefix ${sample}_pb_manta \
        --tmp-dir /scratch/vpagano/tmp
    """

}

process pb_strelka_normalOnly {
    input:
        tuple val(sample), path(bam)
        path reference

    output:
        tuple val("${sample}"), path("${sample}*")
        path "${sample}_pb_strelka.strelka_work"

    when:
        params.pb_strelka

    cpus params.strelkaCpus

    script:
    """
        source /etc/profile.d/modules.sh
        module load parabricks/${params.pb_ver} 
        pbrun strelka --ref ${reference} \
        --in-bams '${bam}' \
        --num-threads ${task.cpus} \
        --out-prefix ${sample}_pb_strelka \
        --tmp-dir /scratch/vpagano/tmp
    """

}
