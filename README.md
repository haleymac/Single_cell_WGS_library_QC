# Single cell WGS library QC

A one-stop shop for all (most?) of your single cell fastq and bam QCing needs.


#### Directories
**fastq_QC**: contains a snakemake workflow to run fastqc, multiqc and fastqscreen to check for contamination in your fastq files, adapter or otherwise

**WGS_BAM_QC**: contains a snakemake workflow to run a number of single cell whole genome BAM library metrics to see if your DLP run worked or not. Most tools run are from GATK, though some are custom 

**scripts**: contains a number of python scripts with helpful functions to run the above workflow, or plot their output


... better explanations to come on this README when my thesis writing is complete...
