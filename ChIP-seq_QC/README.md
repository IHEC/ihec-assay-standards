# Folder contents:

This folder contains a .PDF file with the definitions and a bash script to calculate the ChIP-seq QC metrics.

Note: v0.3 has been tested only on single-end BP ChIP-seq datasets. 

Here's an example of how to run the script:

bash calculate_ChIP-seq_QC_metrics.bash 
     -c $ChIP_original_BAM_file     
     -t H3K4me3 
     -n $ChIP_sampleName 
     -i $Input_sampleName 
     -p $bed_file  
     > ${ChIP_sampleName}_qc_metrics.tsv
