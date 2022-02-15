# TARDIS Documentation

## Initial setup

The TARDIS repository is on github. Simply clone that repo and you can run TARDIS.

    git clone git@github.com:v-pagano/tardis.git

## Running a simple pipeline

Running a pipeline is done through nextflow.  It will also utilize Tower, giving you a nice GUI for run stats and progess. To run a simple pipeline you just run from the command line:

    nextflow run main.nf --input <Your json file>

This will run an alignment-only pipeline using parabricks.

### Using csv or tsv for your input file

TARDIS also supports csv or tsv files for input.

The order for the columns is:

subject
sex
status
sample
lane
fastq1
fastq2
rg

This format is also compatible with nf-core pipelines.

## Configuring for a different cluster or setup

Everything in TARDIS is configurable. Parameters in the nextflow.config file point to all of the containers that are used by TARDIS. If any parameter needs to be tuned for your specific pipeline (# of CPUs for a step, slurm partition to run on, location of reference files), you can create your own copy of the nextflow.config file.  Then run it with this command:

    nextflow -C my_new_config.config run main.nf --input <Your json file>

You can also add a separate profile to the nextflow.config file. There are currently four profiles: dback, gemini, local and coh_apps. You can run TARDIS with a profile like this:

    nextflow run main.nf --input <Your json file> -profile local

You can also run multiple profiles together.  For example to run with the coh_apps containers on dback:

    nextflow run main.nf --input <Your json file> -profile coh_apps,dback

## Command line parameters

All steps in TARDIS are optional and can be put in or removed from the command line.  Most are done with true/false.  Here is a breakdown by pipeline_type.

### Pipeline_type

There are currently 3 different pipeline_types: Alignment, Constitutional and Somatic. You set these with the pipeline_type parameter.

    nextflow run main.nf --pipeline_type Constitutional --input tumor_100K.json

### Mapping_type

There are currently 3 different mapping_types: bwa, bwa2 and parabricks. You set these with the mapping_type parameter.

    nextflow run main.nf --mapping_type bwa2 --input tumor_100K.json

Note that parabricks requires a GPU node and a parabricks license.

### Mapping parameters

Base recalibration is off by default but can be enabled with baserelibration flag.

    nextflow run main.nf --mapping_type parabricks --baserecalibration --input tumor_100K.json

### Variant Callers

There are currently 4 variant callers available: haplotypecaller, deepvariant, pb_haplotypecaller and pb_deepvariant. These are all off by default, but can be turned by adding them to the commmand line.

    nextflow run main.nf --pipeline_type Constitutional --haplotypecaller --pb_deepvariant --input tumor_100K.json

### Annotation

Currently VEP and snpeff are available. These are all off by default, but can be turned by setting them to true on the commmand line.

    nextflow run main.nf --pipeline_type Constitutional --haplotypecaller --vep --input tumor_100K.json

### Using a different reference

TARDIS defaults to a human reference for mapping. To use a different reference change either the reference or pb_reference in the options. These are normally hidden, but you can toggle to make them viewable.

## Publishing results

If you want to copy your output files (BAMs, VCFs etc) to a directory, set the publishResults flag and set the outputFolder.

    nextflow run main.nf --outputFolder /scratch/vpagano/results --publishResults

