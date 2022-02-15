#!/bin/bash -ue
source /etc/profile.d/modules.sh
module load parabricks/3.7.0-1.beta1_V100 
pbrun germline --bwa-options '-K 100000000 -Y' --ref /home/tgenref/canis_familiaris/canfam3.1/canfam3.1_tgen/genome_reference/Canis_familiaris.CanFam3.1.dna.toplevel.fa          --in-fq /scratch/vpagano/play/SRR10351554.sra_1_100K.fastq.gz /scratch/vpagano/play/SRR10351554.sra_2_100K.fastq.gz "@RG\tID:SUPERFQS01_1_K18088\tLB:K18088\tPU:SUPERFQS01_1\tSM:GIAB_NA12878_1_CL_Whole_C1\tPL:ILLUMINA\tCN:tgen\tPM:HiSeq2500\tBC:Unknown"          --out-bam 'SRR10351554.bam'         --out-variants 'SRR10351554_pb_haplotypecaller.vcf'         --tmp-dir /scratch/vpagano/tmp
