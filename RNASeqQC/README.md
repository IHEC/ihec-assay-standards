This folder contains a .PDF document with the definitions and shell script for computation of the RNA-seq QC metrics.

To calculate the metrics please use:

- samtools v 1.3.1

- picard v 2.9.0


## Prerequisites:
For the script to run you have to define the following environment variables:


     export PICARD_290=<path-to-picard-tools-2.9.0>

     export SAMTOOLS_131=<path-to-samtools-1.3.1>

     export $TMP_DIR=<path-to-tmp-dir>

The script will create a folder in the chosen working directory with the name of the bam file. A statistics output file will be written inside in text format. This allows running many samples in parallel in your computational environment. Temporary files (e.g. duplicate marked BAM files) are removed to minimize disk usage.


## How to use:
To run the statistics on your BAM file mapped to the human genome (assembly GRCh37 or GRCh38) you need to provide the script with the BAM file to be analyzed and two BED files (provided in the above subdirectories), one containing intergenic regions and one with rRNA coordinates as the following:


	[sullrich@cluster]$ ./rna_stats.sh

	   USAGE: ./rna_stats.sh <SAMPLE_PATH> <INTERGENIC_BED> <RRNA_BED>
	   The script computes RNA-seq quality metrics according to IHEC standards

	   Input files in the following order are required:
	   1. path to BAM file to be analyzed
	   2. Bed file with intergenic coordinates
	   3. Bed file with rRNA coordinates


## Example:

./rna_stats.sh SAMPLE_ID.bam intergeneic_regions_gencode22_GRCh38.bed rRNA_Mt-rRNA_gencode22_GRCh38.bed

The metrics are found in the text file SAMPLE_ID_read_stats.txt created in the directory SAMPLE_ID and have the following format:

	FRACTION_MAPPED		0.9652343
	FRACTION_INTERGENIC	0.0222358
	FRACTION_RRNA		0.0009536
	FRACTION_DUPLICATED	0.7463131
