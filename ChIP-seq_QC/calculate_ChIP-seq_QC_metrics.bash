#!/bin/bash

usage() { echo "Usage: $0 [-c <ChIP_original_BAM_file>] [-t <H3K27ac|H3K27me3|H3K36me3|H3K4me1|H3K4me3|H3K9me3|Input|H2AFZ|H3ac|H3K4me2|H3K9ac>] [--cn <ChIP_sampleName>] [--in <Input_sampleName>] [-p <bed_file>]" 1>&2; exit 1; }

while getopts ":c:t:n:i:p:" o; do
    case "${o}" in
        c)
            c=${OPTARG}
            ;;
        t)
            t=${OPTARG}
            ((t == "H3K27ac" || t == "H3K27me3" || t == "H3K36me3" || t == "H3K4me1" || t == "H3K4me3" || t == "H3K9me3" || t == "Input" || t == "H2AFZ" || t == "H3ac" || t == "H3K4me2" || t == "H3K9ac" )) || usage
            ;;
        n)
            cname=${OPTARG}
            ;;
        i)
            iname=${OPTARG}
            ;;
        p)
            p=${OPTARG}
            ;;
        *)
            usage
            ;;
    esac
done
shift $((OPTIND-1))

if [ -z "${c}" ] || [ -z "${t}" ] || [ -z "${cname}" ] || [ -z "${iname}" ] || [ -z "${p}" ]; then
    usage
fi

export PATH=${PATH}:/homes/kostadim/bin/:/nfs/1000g-work/G1K/work/bin/bedtools2-2.20.1/bin/

## The following commands assume that there is a pair of BAM files, one for the ChIP and one for the Input, $ChIP_original_BAM_file and $Input_original_BAM_file for two samples labelled, $ChIP_sampleName and $Input_sampleName, respectively.
## The following steps are shown for the $ChIP_sampleName but have to be applied to the $Input_sampleName too:

## Sort the BAM file by coordinate
if [[ ! -s ${cname}_original.sorted.bam ]]
then
  java -Xmx2048m -jar /nfs/1000g-work/G1K/work/bin/picard-tools-1.137/picard.jar SortSam INPUT=$c OUTPUT=${cname}_original.sorted.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT
fi

## Mark, but not remove, duplicate reads
if [[ ! -s ${cname}_markDup.bam ]]
then
  java -Xmx2048m -jar /nfs/1000g-work/G1K/work/bin/picard-tools-1.137/picard.jar MarkDuplicates INPUT=${cname}_original.sorted.bam OUTPUT=${cname}_markDup.bam METRICS_FILE=${cname}_original.sorted_metrics.out REMOVE_DUPLICATES=false ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT
fi

## Remove unmapped reads and those with mapping quality less than 5:
if [[ ! -s ${cname}_quality_filtered.bam ]]
then
  samtools view -b -F 4 -q 5 ${cname}_markDup.bam > ${cname}_quality_filtered.bam
fi

## Remove duplicate reads and those with mapping quality less than 5:
if [[ ! -s ${cname}_dedup.bam ]]
then 
  ## Remove duplicate reads:
  samtools view -b -F 1024 ${cname}_quality_filtered.bam > ${cname}_dedup.bam

  ## Index the final deduplicated BAM file
  samtools index ${cname}_dedup.bam
fi

#2.	Mappability
#We want to extract three mapping statistics:
## The original number of reads and the number of those aligned:
if [[ ! -s ${cname}_original_flagstat.txt ]]
then
  samtools flagstat ${cname}_markDup.bam > ${cname}_markDup_flagstat.txt
fi

total_reads=`grep "in total" ${cname}_markDup_flagstat.txt | sed -e 's/ + [[:digit:]]* in total .*//'`
mapped_reads=`grep "mapped (" ${cname}_markDup_flagstat.txt | sed -e 's/ + [[:digit:]]* mapped (.*)//'`
dupped_reads=`grep "duplicates" ${cname}_markDup_flagstat.txt | sed -e 's/ + [[:digit:]]* duplicates$//'`
dup_rate=$(echo "${dupped_reads}/${mapped_reads}" | bc -l)

## Finally, the number of singletons for paired-end data sets can be calculated using:
left_singletons=`grep "singletons" ${cname}_markDup_flagstat.txt | sed -e 's/ + [[:digit:]]* singletons .*//'`
right_singletons=`grep "singletons" ${cname}_markDup_flagstat.txt | sed -e 's/[[:digit:]]* + //;s/ singletons .*//'`
singletons=$((left_singletons+right_singletons))

## The final number of reads:
if [[ ! -s ${cname}_dedup_flagstat.txt ]]
then
    samtools flagstat ${cname}_dedup.bam > ${cname}_dedup_flagstat.txt
fi

final_reads=`grep "mapped (" ${cname}_dedup_flagstat.txt | sed -e 's/ + [[:digit:]]* mapped (.*)//'`

#3.	Calculating Jensen-Shannon distance (JSD)

#To calculate the Jensen-Shannon distance we run:
## Attention: Regarding the bin size (specified in the command below by the ‘-bs’ option) there hasn’t been an agreement on what the optimal bin size is yet. There have been discussions on adopting smaller bin sizes for the sharp peaks and larger bin sizes for the broad peaks.
## No need to remove the blacklisted regions for the JSD calculation.

if [[ ! -s ${cname}_fingerprint.txt ]]
then
  if [[ type == "H3K27ac" || type == "H3K4me3" || type == "H2AFZ" || type == "H3ac" || type == "H3K4me2" || type == "H3K9ac" ]]
  then
    bin_size=200
  else
    bin_size=1000
  fi

  plotFingerprint -b ${cname}_dedup.bam ${iname}_dedup.bam -bs ${bin_size} -l ${cname} ${iname} --JSDsample ${iname}_dedup.bam --outQualityMetrics ${cname}_fingerprint.txt -plot ${cname}_fingerprint.png -p 8

fi

js_dist=`grep ${cname} ${cname}_fingerprint.txt | cut -f 8`
chance_div=`grep ${cname} ${cname}_fingerprint.txt | cut -f 12` 


#4.     Calculating FRiP scores

reads_under_peaks=`bedtools intersect -wa -bed -abam ${cname}_dedup.bam -b ${p} | wc -l`
frip=$(echo "${reads_under_peaks}/${final_reads}" | bc -l)

printf "ChIP_name\tInput_name\ttotal_reads\tmapped_reads\tdupped_reads\tdup_rate\tsingletons\tfinal_reads\tjs_dist\tchance_div\tfrip\n" 
printf "%s\t%s\t%d\t%d\t%d\t%.4f\t%d\t%d\t%.4f\t%.4f\t%.4f\n" "$cname" "$iname" "$total_reads" "$mapped_reads" "$dupped_reads" "$dup_rate" "$singletons" "$final_reads" "$js_dist" "$chance_div" "$frip"


