#!/bin/bash

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

SAMPLE_PATH=$1
INTERGENIC_BED=$2
RRNA_BED=$3

SAMPLE_NAME=$(basename $SAMPLE_PATH .bam)

mkdir $SAMPLE_NAME

cd $SAMPLE_NAME

#computation of IHEC RNA-seq stats
#get the fraction of mapped reads
MAPPED_READS=$(samtools view -F 0x904 -c $SAMPLE_PATH)

samtools view -F 0x900 -c $SAMPLE_PATH | awk -v MAPPED_READS=$MAPPED_READS '{print "FRACTION_MAPPED\t"MAPPED_READS/$1}' > ${SAMPLE_NAME}_read_stats.txt

#get the number of reads falling into intergenic regions with samtools
samtools view -c -F 0x900 $SAMPLE_PATH -L $INTERGENIC_BED | awk -v MAPPED_READS=$MAPPED_READS '{print "FRACTION_INTERGENIC\t"$1/MAPPED_READS}' >> ${SAMPLE_NAME}_read_stats.txt

#get the number of reads from rRNA with samtools
samtools view -c -F 0x900 $SAMPLE_PATH -L $RRNA_BED | awk -v MAPPED_READS=$MAPPED_READS '{print "FRACTION_RRNA\t"$1/MAPPED_READS}' >> ${SAMPLE_NAME}_read_stats.txt

#get number of duplicates
java -Xmx4g -jar $PICARD/picard.jar MarkDuplicates VALIDATION_STRINGENCY=LENIENT I=${SAMPLE_NAME}_no_multimap.bam O=${SAMPLE_NAME}_noMULTI_noDUP.bam M=${SAMPLE_NAME}_duplicated.txt

awk -F $'\t' '/PERCENT_DUPLICATION/{getline; print "FRACTION_DUPLICATED\t"$8}' ${SAMPLE_NAME}_duplicated.txt >> ${SAMPLE_NAME}_read_stats.txt

rm ${SAMPLE_NAME}_noMULTI_noDUP.bam
rm ${SAMPLE_NAME}_duplicated.txt
