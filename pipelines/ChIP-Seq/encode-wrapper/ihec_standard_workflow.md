# IHEC ChIP-Seq standard workdlows

See the ENCODE reference for input format: https://github.com/ENCODE-DCC/chip-seq-pipeline2/blob/master/docs/input.md

* For 75bp PET (histone) ChIP-Seq:

See the configuration see example `cemt0007_h3k4me3_mnt_ext_0.json` provided. 

* For readlengths longer than 75bp PET (histone) ChIP-Seq:

Trim fastq to 75bp pre-alignment using trimfastq.py from the singularity pipeline image (also: https://github.com/ENCODE-DCC/chip-seq-pipeline2/blob/master/src/trimfastq.py). Use 

    Singularity chip_seq_pipeline_v1_1_2.simg:~> python /software/chip-seq-pipeline/src/trimfastq.py 
        usage: python /software/chip-seq-pipeline/src/trimfastq.py <inputfilename> <bpToKeep | max> [-trim5 bp] [-flowcellID flowcell] [-addEnd 1 | 2] [-replace string newstring | blank] [-renameIDs prefix] [-stdout]
	      the -trim5 option will trim additional bp from the 5 end, i.e. if you want the middle 36bp of 38bp reads, use 36 as bp to keep and 1 as the trim5 argument
	      Use - to specify standard input, the script will print(to standard output by default
	      The script can read compressed files as long as they have the correct suffix - .bz2 or .gz

* For SET, just set `"chip.paired_end" : false.` 
