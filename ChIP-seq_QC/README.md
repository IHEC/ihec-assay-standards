# Folder contents:

This folder contains a .PDF file with the definitions and a bash script to calculate the ChIP-seq QC metrics.

Note: v0.3 has only been tested on single-end ChIP-seq datasets. 

Here's an example of how to run the script:

bash calculate_ChIP-seq_QC_metrics.bash

     -c $ChIP_original_BAM_file     

     -t H3K4me3 
     
     -n $ChIP_sampleName 
     
     -i $Input_sampleName 
     
     -p $bed_file  
     
     > ${ChIP_sampleName}_qc_metrics.tsv

To calculate the metrics please use:

- samtools v 1.3.1

- picard v 2.9.0

- plotFingerprint v 2.4.2 (deepTools)

- bedtools v 2.26
