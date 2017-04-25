This folder contains a .PDF file with the definitions and two bash scripts: one to pre-process the BAM files and one to calculate the ChIP-seq QC metrics.

To calculate the metrics please use:

- samtools v 1.3.1

- picard v 2.9.0

- plotFingerprint v 2.5.0.1 (deepTools)


## Prerequisites:
For the scripts to run you have to define the following environment variables:


     export PICARD_290=<path-to-picard-tools-2.9.0>

     export DEEPTOOLS_2501=<path-to-deeptools-2.5.0.1>

     export SAMTOOLS_131=<path-to-samtools-1.3.1>


Optionally, you should also set the working directory, using:


     export $WORKING_DIR=<path-to-working-directory>


If $WORKING_DIR is not set, then the scripts will write any output to the current working directory.


     export $TMP_DIR=<path-to-tmp-dir>


If $TMP_DIR is set, then a dataset-specific subfolder will be created that will be used when running picard tools.

## How to use:
The first step is to pre-process all ChIP-seq datasets, using `preprocess_BAM_files.bash`.


     bash preprocess_BAM_files.bash \
          
          -c $ChIP_original_BAM_file \

          -f $ChIP_sampleName 


This will create a directory per dataset: $WORKING_DIR/$ChIP_sampleName/ where the following files per dataset will be stored:

+ ${ChIP_sampleName}_original.sorted.bam: if the original BAM file is unsorted, a new sorted BAM file will be created
+ ${ChIP_sampleName}_markDup.bam: a BAM file where duplicates are marked
+ ${ChIP_sampleName}_original.sorted_metrics.out: the METRICS_FILE reported by MarkDuplicates
+ ${ChIP_sampleName}_quality_filtered.bam: the final BAM file without duplicates, unampped reads  or reads with mapping quality < 5
+ ${ChIP_sampleName}_dedup.bam.bai

The second step is to calculate the ChIP-seq quality statistics, using `calculate_ChIP-seq_QC_metrics.bash`

     bash calculate_ChIP-seq_QC_metrics.bash \

          -c $ChIP_original_BAM_file \

          -t <H3K27ac|H3K27me3|H3K36me3|H3K4me1|H3K4me3|H3K9me3|Input|H2AFZ|H3ac|H3K4me2|H3K9ac> \
     
          -f $ChIP_sampleName \ 
     
          -u $Input_sampleName \ 
     
          -p $bed_file \ ## Optional, only provide it when the  set is NOT an Input

          -n no_threads \ ## Optional, by default is set to 1
     
          > ${ChIP_sampleName}_qc_metrics.tsv


Attention: This script will be looking for the files generated using the `preprocess_BAM_files.bash`, so make sure to not chnage the $WORKING_DIR between scripts. 

## Output:

The ChIP-seq metrics are reported in a text file `$ChIP_sampleName_read_stats.txt` created in the directory `$ChIP_sampleName`.

## Example:

Suppose we have two ChIP-seq datasets to process:
1. a H3K4me3 with the BWA output being stored in a file called: `S00XDKH1.ERX712764.H3K4me3.unfiltered.bwa.GRCh38.20150503.bam` and the peaks stored in a file called: `S00XDKH1.ERX712764.H3K4me3.bwa.GRCh38.20150527.bed.gz`
2. an Inout with the BWA outut being stored in a file caled: `S00XDKH1.ERX941048.Input.unfiltered.bwa.GRCh38.20150503.bam` 

We would first preprocess both datasets using the two following commands:


     bash preprocess_BAM_files.bash \
          -c S00XDKH1.ERX712764.H3K4me3.unfiltered.bwa.GRCh38.20150503.bam \
          -f S00XDKH1.ERX712764.H3K4me3

     bash preprocess_BAM_files.bash \
          -c S00XDKH1.ERX941048.Input.unfiltered.bwa.GRCh38.20150503.bam \
          -f S00XDKH1.ERX941048.Input
     

And then calculate the ChIP-seq statistics using the following commands:


     bash calculate_ChIP-seq_QC_metrics.bash \
          -t H3K4me3 -n 8 \
          -f S00XDKH1.ERX712764.H3K4me3 \
          -u S00XDKH1.ERX941048.Input \
          -p S00XDKH1.ERX712764.H3K4me3.bwa.GRCh38.20150527.bed.gz


     bash calculate_ChIP-seq_QC_metrics.bash \
          -t Input \
          -f S00XDKH1.ERX941048.Input \
          -u S00XDKH1.ERX941048.Input



