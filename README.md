# TARDIS Documentation

## Initial setup

The TARDIS repository is on github. Simply clone that repo and you can run TARDIS.

    git clone git@github.com:v-pagano/tardis.git

## Running a simple pipeline

First thing, you need to let the pipeline know what your input format is.  There are currently 4 possible inputs: fastq, ubam, bam or vcf.  You set this with the inputType parameter like so:

    nextflow run main.nf --inputType fastq

Once that is set, simply point to your input files:

    nextflow run main.nf --inputType bam --input "/scratch/vpagano/play/*.bam"

Paired fastqs require you to tell nextflow where the _1 or _2 is located in the filename, you do this with the {1,2} setup:

    nextflow run main.nf --inputType fastq --input "/scratch/vpagano/play/*L1R{1,2}*"

Once this done, you simply have to decide what you want to run in the pipeline, for example to run parabricks fq2bam on some fastq files:

    nextflow run main.nf --inputType fastq --input "/scratch/vpagano/play/*L1R{1,2}*" --parabricks

## Pipeline parameters

### Pipeline steps

You control which steps that the pipeline should run with command line parameters, they are simply on/off toggles.  For example to run parabricks haplotypecaller on a bam file:

    nextflow run main.nf --inputType bam --input "/scratch/vpagano/play/*.bam" --pb_haplotypecaller

You can string together as many variations as you would like. For example to run bwamem2 alignment with parabricks deepvariant, VEP annotation and fastqc:

    nextflow run main.nf --inputType fastq --input "/scratch/vpagano/play/*L1R{1,2}*" --bwa2 --pb_deepvariant --vep --fastqc

#### Fastq Prep & QC

| Step                               | Parameter           |
|------------------------------------|---------------------|
| Run fastqc                         | --fastqc            |
| Create a uBAM from the fastq files | --ubam              |
| Run parabricks multiQC             | --pb_collectmetrics |

#### Mapping and alignment

| Step               | Parameter    |
|--------------------|--------------|
| Parabricks fq2bam  | --parabricks |
| BWA MEM alignment  | --bwa        |
| BWA MEM2 alignment | --bwa2       |

#### Variant Callers

| Step                       | Parameter            |
|----------------------------|----------------------|
| Parabricks haplotypecaller | --pb_haplotypecaller |
| Parabricks deepvariant     | --pb_deepvariant     |

#### Annotation

| Step   | Parameter |
|--------|-----------|
| SNPEff | --snpeff  |
| VEP    | --vep     |

#### File output

| Step                         | Parameter         |
|------------------------------|-------------------|
| Publish all files            | --publishResults  |
| Where to output the files    | --outputFolder    |

#### Other Settings

| Step                    | Parameter           |
|-------------------------|---------------------|
| Produce gvcf            | --gvcf              |
| Use BQSR during mapping | --baserecalibration |

### Profiles

Profiles allow you to change many command line parameters with a simple profile. For example to run with the local profile

    nextflow run main.nf -profile local

Here are the current profiles:

| Profile | Modifications                                 |
|---------|-----------------------------------------------|
| dback   | Configured for the dback cluster              |
| gemini  | Configured for the gemini cluster             |
| local   | Configured to run on a single local node      |
| k9      | Configured to use boxer Tasha reference       |
| reports | Outputs reports on CPU usage and elapsed time |
