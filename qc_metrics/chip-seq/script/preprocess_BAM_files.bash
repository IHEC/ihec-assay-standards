#!/bin/bash

set -euf -o pipefail

usage() {
  echo "Usage: "`basename $0`" [-c <ChIP_original_BAM_file> ] "
  echo "          [-f <new_files_prefix>]"
 }

while getopts "c:t:f:i::u:p::n::" o; do
    case "${o}" in
        c)
            c=${OPTARG}
            ;;
        f)
            cname=${OPTARG}
            ;;
        *)
            ;;
    esac
done

if [[ -z "${c+x}" || -z "${cname+x}" ]]
then
  usage
  exit 1
fi

if [[ -z ${PICARD_290+x} || -z ${SAMTOOLS_131+x} || -z ${DEEPTOOLS_2501+x} ]]; then
  echo "Environment variable PICARD_290 (path to picard jar v2.9.0), " >&2
  echo "SAMTOOLS_131 (path to samtools v1.3.1) and " >&2
  echo "DEEPTOOLS_2501 (path to deeptools v2.5.0.1) must be set" >&2
  exit 1
fi

## Set the working directory to current directory, if one hasn't been set.
if [[ -z ${WORKING_DIR+x} ]]; then
  WORKING_DIR="."
  echo "WORKING_DIR unset; used default: $WORKING_DIR"
fi
mkdir -p $WORKING_DIR/$cname || true
cd $WORKING_DIR/$cname

## Create dataset-speciic TMP_DIR for picard, if $TMP_DIR is set.
if [[ -z $TMP_DIR ]]; then
  PICARD_TMP_DIR=""
else
  TEMP="$TMP_DIR/${cname}/temp"
  mkdir -p $TEMP || true
  PICARD_TMP_DIR="TMP_DIR=$TEMP"
fi  


sortBamIfNeeded(){
  bam=$1
  sortedBam=$2
  set +e
  $SAMTOOLS_131/samtools view -H $bam | grep 'SO:coordinate'
  IS_SORTED=$?
  echo "isSorted:$IS_SORTED $bam"
  set -e
  if [[ $IS_SORTED -ne 0 ]]; then
    echo "sorting.."
    java -Xmx2048m -jar $PICARD_290/picard.jar SortSam INPUT=$bam OUTPUT=$sortedBam SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT $PICARD_TMP_DIR
    echo "return[sorting]:$? "
    file_to_mark_duplicates=$sortedBam
  else
    echo "$bam sorted"
    file_to_mark_duplicates=$bam
  fi
}


sortBamIfNeeded $c ${cname}_original.sorted.bam
echo "duplicateMarking:$file_to_mark_duplicates"


## Mark, but not remove, duplicate reads
if [[ ! -s ${cname}_markDup.bam ]]
then
  java -Xmx2048m -jar $PICARD_290/picard.jar MarkDuplicates INPUT=$file_to_mark_duplicates OUTPUT=${cname}_markDup.bam METRICS_FILE=${cname}_original.sorted_metrics.out REMOVE_DUPLICATES=false ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT $PICARD_TMP_DIR
  echo "return[marked duplicates]:$? "
else
  echo "duplicate marked bam exists: ${cname}_markDup.bam"
fi

## Remove unmapped read, duplicate reads and those with mapping quality less than 5:
${SAMTOOLS_131}/samtools view -b -F 3844 -q 5 ${cname}_markDup.bam > ${cname}_dedup.bam
echo "return[F3844q5 filtered]:$? "

## Index the final deduplicated BAM file
${SAMTOOLS_131}/samtools index ${cname}_dedup.bam
echo "return[indexed]:$? ${cname}_dedup.bam"

