#!/bin/bash

set -euf -o pipefail 

if [ $# -lt 3 ]; then
  echo
  echo "   USAGE: $0 <SAMPLE_PATH> <INTERGENIC_BED> <RRNA_BED> [<THREADS>]"
  echo "   The script computes RNA-seq qualtity metrics according to IHEC standards"
  echo ""
  echo "   Input files in the following order are required:"
  echo "   1. path to BAM file to be analyzed"
  echo "   2. Bed file with intergenic coordinates"
  echo "   3. Bed file with rRNA coordinates"
  echo
  echo "   An optional 4th argument can be used to specify the maximum number of threads to use (default: 1)"
  echo
  exit 1
fi

if [[ -z ${SAMBAMBA+x} ]]; then 
  echo "Environment variable SAMBAMBA (path to sambamba binary) must be set"
  exit 1
fi

resolve() {
    if [[ "$1" = /* ]];then
        echo $1
        return
    fi
    HEAD=$1
    if [ ! -d $1 ]; then
        HEAD=$(dirname $1)
        TAIL="/$(basename $1)"
    fi
    cd $HEAD && echo "$(pwd -P)${TAIL-}"
}

SAMPLE_PATH=$(resolve $1)
INTERGENIC_BED=$(resolve $2)
RRNA_BED=$(resolve $3)
THREADS=${4-"1"}

SAMPLE_NAME=$(basename $SAMPLE_PATH .bam)

WORKING_DIR=${QC_DIR:-"."} # by default working directory is current
mkdir -p $WORKING_DIR/$SAMPLE_NAME || true
cd $WORKING_DIR/$SAMPLE_NAME

BASE_FILTER='not (secondary_alignment or supplementary)'
MAPPED_FILTER='not (secondary_alignment or supplementary or unmapped)'

# computation of IHEC RNA-seq stats
# get the fraction of mapped reads
MAPPED_READS=$($SAMBAMBA view -c -t $THREADS -F "$MAPPED_FILTER" $SAMPLE_PATH) 

N_NO_MULTIMAP=$($SAMBAMBA view -c -t $THREADS -F "$BASE_FILTER" $SAMPLE_PATH)

echo $N_NO_MULTIMAP | awk -v MAPPED_READS=$MAPPED_READS '{print "FRACTION_MAPPED\t"MAPPED_READS/$1}' > ${SAMPLE_NAME}_read_stats.txt

# get the number of reads falling into intergenic regions
$SAMBAMBA view -t $THREADS -c -F "$BASE_FILTER" $SAMPLE_PATH -L $INTERGENIC_BED | awk -v MAPPED_READS=$MAPPED_READS '{print "FRACTION_INTERGENIC\t"$1/MAPPED_READS}' >> ${SAMPLE_NAME}_read_stats.txt

# get the number of reads from rRNA regions
$SAMBAMBA view -t $THREADS -c -F "$BASE_FILTER" $SAMPLE_PATH -L $RRNA_BED | awk -v MAPPED_READS=$MAPPED_READS '{print "FRACTION_RRNA\t"$1/MAPPED_READS}' >> ${SAMPLE_NAME}_read_stats.txt

# get number of duplicates
# add TMPDIR to mark dups
TEMP="/tmp"
if [[ ! -z ${TMPDIR+x} ]]; then
  TEMP="$TMPDIR/$SAMPLE_NAME/temp"
  mkdir -p $TEMP || true
fi  
$SAMBAMBA markdup -t $THREADS $SAMPLE_PATH --tmpdir ${TEMP} ${SAMPLE_NAME}_noMULTI_noDUP.bam
$SAMBAMBA view -t 8 -c -F 'duplicate' ${SAMPLE_NAME}_noMULTI_noDUP.bam | awk -v MAPPED_READS=$MAPPED_READS '{print "FRACTION_DUPLICATED\t"$1/MAPPED_READS}' >> ${SAMPLE_NAME}_read_stats.txt
# cleanup
rm ${SAMPLE_NAME}_noMULTI_noDUP.bam{,.bai}
