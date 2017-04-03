#!/bin/bash

usage() {
  echo "Usage: "`basename $0`" [-c <ChIP_original_BAM_file> ] "
  echo "          [-f <new_files_prefix>]"
  exit 1
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

if [[ -z "${c}" || -z "${cname}" ]]
then
    usage
fi

if [[ -z ${PICARD_290+x} || -z ${SAMTOOLS_131+x} || -z ${DEEPTOOLS_242+x} ]]; then
  echo "Environment variable PICARD_290 (path to picard jar v2.9.0), " >&2
  echo "SAMTOOLS_131 (path to samtools v1.3.1) and " >&2
  echo "DEEPTOOLS_242 (path to deeptools v2.4.2) must be set" >&2
  exit 1
fi

## Set the working directory to current directory, if one hasn't been set.
if [[ -z ${WORKING_DIR} ]];
then
  echo "The WORKING_DIR variable isn't set. The current directory will be used instead." >&2 
  WORKING_DIR=${QC_DIR:-"."}
fi
mkdir -p $WORKING_DIR/$cname || true
cd $WORKING_DIR/$cname

## Sort the original BAM file by coordinate, if it isn't yet.

if ! $SAMTOOLS_131/samtools view -H $file2process | grep -q "SO:coordinate"
then
  java -Xmx2048m -jar $PICARD_290/picard.jar SortSam INPUT=$c OUTPUT=${cname}_original.sorted.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT

  file_to_mark_duplicates=${cname}_original.sorted.bam
else
  file_to_mark_duplicates=$c
fi

## Mark, but not remove, duplicate reads
if [[ ! -s ${cname}_markDup.bam ]]
then
  java -Xmx2048m -jar $PICARD_290/picard.jar MarkDuplicates INPUT=$file_to_mark_duplicates OUTPUT=${cname}_markDup.bam METRICS_FILE=${cname}_original.sorted_metrics.out REMOVE_DUPLICATES=false ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT
fi

## Remove unmapped read, duplicate reads and those with mapping quality less than 5:
${SAMTOOLS_131}/samtools view -b -F 1028 -q 5 ${cname}_markDup.bam > ${cname}_dedup.bam

## Index the final deduplicated BAM file
${SAMTOOLS_131}/samtools index ${cname}_dedup.bam


