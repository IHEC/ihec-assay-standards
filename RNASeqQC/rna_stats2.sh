#!/bin/bash

set -eufx -o pipefail # -x is noisy but keep on till in production 

if [ $# -lt 3 ]; then
  echo
  echo "   USAGE: $0 <SAMPLE_PATH> <INTERGENIC_BED> <RRNA_BED>"
  echo "   The script computes RNA-seq qualtity metrics according to IHEC standards"
  echo ""
  echo "   Input files in the following order are required:"
  echo "   1. path to BAM file to be analyzed"
  echo "   2. Bed file with intergenic coordinates"
  echo "   3. Bed file with rRNA coordinates"
  echo
  exit 1
fi

if [[ -z ${PICARD_290+x} || -z ${SAMTOOLS_131+x} ]]; then 
  echo "Environment variable PICARD_290 (path to picard jar v2.9.0) and SAMTOOLS_131 (path to samtools v1.3.1) must be set"
  exit 1
fi

SAMPLE_PATH=$1
INTERGENIC_BED=$2
RRNA_BED=$3

SAMPLE_NAME=$(basename $SAMPLE_PATH .bam)

WORKING_DIR=${QC_DIR:-"."} # by default working directory is current
mkdir -p $WORKING_DIR/$SAMPLE_NAME || true
cd $WORKING_DIR/$SAMPLE_NAME

#computation of IHEC RNA-seq stats
#get the fraction of mapped reads
$SAMTOOLS_131/samtools view -F 0x904 -c $SAMPLE_PATH > $SAMPLE_NAME.n_mapped #&
#PID=$! # send this process into background and wait on it before reading $SAMPLE_NAME.n_mapped

# create 0x900 filtered bam and count those reads at the same time with tee redireting to no_multimap_bam
# faster because sone concurrently + no uncompressed data streamed
N_NO_MULTIMAP=$($SAMTOOLS_131/samtools view -b -F 0x900 $SAMPLE_PATH | tee ${SAMPLE_NAME}_no_multimap.bam | $SAMTOOLS_131/samtools view -c -)

# wait on $SAMPLE_NAME.n_mapped process
#wait $PID 
MAPPED_READS=$(cat $SAMPLE_NAME.n_mapped) 
echo $N_NO_MULTIMAP | awk -v MAPPED_READS=$MAPPED_READS '{print "FRACTION_MAPPED\t"MAPPED_READS/$1}' > ${SAMPLE_NAME}_read_stats.txt

#get the number of reads falling into intergenic regions with samtools
$SAMTOOLS_131/samtools view -c ${SAMPLE_NAME}_no_multimap.bam -L $INTERGENIC_BED | awk -v MAPPED_READS=$MAPPED_READS '{print "FRACTION_INTERGENIC\t"$1/MAPPED_READS}' >> ${SAMPLE_NAME}_read_stats.txt

#get the number of reads from rRNA with samtools
$SAMTOOLS_131/samtools view -c ${SAMPLE_NAME}_no_multimap.bam -L $RRNA_BED | awk -v MAPPED_READS=$MAPPED_READS '{print "FRACTION_RRNA\t"$1/MAPPED_READS}' >> ${SAMPLE_NAME}_read_stats.txt

#get number of duplicates
PICARD_MARK_DUP_CMD="java -Xmx4g -jar $PICARD_290/picard.jar MarkDuplicates VALIDATION_STRINGENCY=LENIENT I=${SAMPLE_NAME}_no_multimap.bam O=${SAMPLE_NAME}_noMULTI_noDUP.bam M=${SAMPLE_NAME}_duplicated.txt"
#add -TMP_DIR argument to mark dups so avoid running out of space. 
if [[ -z $TMP_DIR ]]; then
  $PICARD_MARK_DUP_CMD
else
  TEMP="$TMP_DIR/$SAMPLE_NAME/temp"
  mkdir -p $TEMP || true
  $PICARD_MARK_DUP_CMD TMP_DIR=$TEMP
fi  
awk -F $'\t' '/PERCENT_DUPLICATION/{getline; print "FRACTION_DUPLICATED\t"$9}' ${SAMPLE_NAME}_duplicated.txt >> ${SAMPLE_NAME}_read_stats.txt

# keep everything while debugging  
#rm ${SAMPLE_NAME}_noMULTI_noDUP.bam 
#rm ${SAMPLE_NAME}_duplicated.txt
