#!/bin/bash

usage() { 
  echo "Usage: "`basename $0`" [-t <H3K27ac|H3K27me3|H3K36me3|H3K4me1|H3K4me3|H3K9me3|Input|H2AFZ|H3ac|H3K4me2|H3K9ac>]"
  echo "          [-f <ChIP_file_prefix>]"
  echo "          [-u <Input_file_prefix]" 
  echo "          [-p <ChIP_bed_file>]" 
  echo "          [-n <threads>]"
  exit 1
 }

## By default the number of threads is set to 1;
n=1

while getopts "t:f:u:p::n::" o; do
    case "${o}" in
        t)
            t=${OPTARG}
            ;;
        f)
            cname=${OPTARG}
            ;;
        u)
            iname=${OPTARG}
            ;;
        p)
            p=${OPTARG}
            ;;
        n)
            n=${OPTARG}
            ;;
        *)
            usage
            ;;
    esac
done

## Set the working directory to current directory, if one hasn't been set.
if [[ -z ${WORKING_DIR+x} ]]; then
  WORKING_DIR="."
  echo "WORKING_DIR unset; used default: $WORKING_DIR"
fi

## If the $WORKING_DIR/Input_prefix/Input_prefix_dedup.bam file doesn't exist, then throw an error 
if [[  ! -s $WORKING_DIR/${iname}/${iname}_dedup.bam ]]
then
  echo "ERROR: The $WORKING_DIR/${iname}/${iname}_dedup.bam file doesn't exist or is empty." >&2
  echo "Make sure that the matched deduplicated Input BAM file is $WORKING_DIR/${iname}/${iname}_dedup.bam. Exiting now.." >&2
  exit 1;
fi

if [[ ( -z "${p}" || ! -s $p ) && ! "$t" == "Input" ]]
then
  echo "ERROR: Your sample isn't of type Input but you're not providing a BED file with the peaks, or the file doesn't exist or is empty." >&2
  echo "Run again providing the correct BED file or change experiment type. Exiting now.." >&2
  exit 1;
fi

if [[ -z "${t}" || -z "${cname}" || -z "${iname}" ]]
then
  usage
fi

if [[ ! ( "$t" == "H3K27ac" || "$t" == "H3K27me3" || "$t" == "H3K36me3" || "$t" == "H3K4me1" || "$t" == "H3K4me3" || "$t" == "H3K9me3" || "$t" == "Input" || "$t" == "H2AFZ" || "$t" == "H3ac" || "$t" == "H3K4me2" || "$t" == "H3K9ac" ) ]]
then
  echo "The experiment type defined isn't one of the following." >&2
  echo "H3K27ac | H3K27me3 | H3K36me3 | H3K4me1 | H3K4me3 | H3K9me3 | Input | H2AFZ | H3ac | H3K4me2 | H3K9ac" >&2
fi 

if [[ -z ${PICARD_290+x} || -z ${SAMTOOLS_131+x} || -z ${DEEPTOOLS_2501+x} ]]; then 
  echo "Environment variable PICARD_290 (path to picard jar v2.9.0), " >&2
  echo "SAMTOOLS_131 (path to samtools v1.3.1) and " >&2
  echo "DEEPTOOLS_2501 (path to deeptools v2.5.0.1) must be set" >&2
  exit 1
fi

if [[ ! -s $WORKING_DIR/$cname/${cname}_dedup.bam ]]
then
  echo "ERROR: File $WORKING_DIR/$cname/${cname}_dedup.bam doesn't exist or is empty. Your sample doesn't seem to have been preprocessed yet." >&2
  echo "Please, run preprocess_BAM_files.bash before running this script. Exiting now.." >&2
  exit 1;
fi

if [[  ! -s $WORKING_DIR/$iname/${iname}_dedup.bam ]]
then
  echo "ERROR: File $WORKING_DIR/$iname/${iname}_dedup.bam doesn't exist or is empty. The matched Input sample you provided doesn't seem to have been preprocessed yet." >&2
  echo "Please, run preprocess_BAM_files.bash on the matched Input sample before running this script. Exiting now.." >&2
  exit 1;
fi

#2.	Mappability
#We want to extract three mapping statistics:
## The original number of reads and the number of those aligned:
if [[ ! -s $WORKING_DIR/${cname}/${cname}_markDup_flagstat.txt ]]
then
  ${SAMTOOLS_131}/samtools flagstat $WORKING_DIR/${cname}/${cname}_markDup.bam > $WORKING_DIR/${cname}/${cname}_markDup_flagstat.txt

  echo "return[run samtools flagstat]:$? "
else
  echo "samtools flagstat .txt exists: ${cname}_markDup_flagstat.txt"
fi

supplementarysecondary_reads=`bc <<< $(grep "secondary" $WORKING_DIR/${cname}/${cname}_markDup_flagstat.txt | sed -e 's/ + [[:digit:]]* secondary.*//')+$(grep "supplementary" $WORKING_DIR/${cname}/${cname}_markDup_flagstat.txt | sed -e 's/ + [[:digit:]]* supplementary.*//')`
total_reads=`bc <<< $(grep "in total" $WORKING_DIR/${cname}/${cname}_markDup_flagstat.txt | sed -e 's/ + [[:digit:]]* in total .*//')-$supplementarysecondary_reads`
mapped_reads=`bc <<< $(grep "mapped (" $WORKING_DIR/${cname}/${cname}_markDup_flagstat.txt | sed -e 's/ + [[:digit:]]* mapped (.*)//')-$supplementarysecondary_reads`
dupped_reads=`grep "duplicates" $WORKING_DIR/${cname}/${cname}_markDup_flagstat.txt | sed -e 's/ + [[:digit:]]* duplicates$//'`
dup_rate=$(echo "${dupped_reads}/${mapped_reads}" | bc -l)

## Finally, the number of singletons for paired-end data sets can be calculated using:
singletons=`grep "singletons" $WORKING_DIR/${cname}/${cname}_markDup_flagstat.txt | sed -e 's/ + [[:digit:]]* singletons .*//'`

final_reads=`${SAMTOOLS_131}/samtools flagstat $WORKING_DIR/${cname}/${cname}_dedup.bam | grep "mapped (" | sed -e 's/ + [[:digit:]]* mapped (.*)//'`

#3.	Calculating Jensen-Shannon distance (JSD)

#To calculate the Jensen-Shannon distance we run:
## Attention: Regarding the bin size (specified in the command below by the ‘-bs’ option) there hasn’t been an agreement on what the optimal bin size is yet. There have been discussions on adopting smaller bin sizes for the sharp peaks and larger bin sizes for the broad peaks.
## No need to remove the blacklisted regions for the JSD calculation.

if [[ ! -s $WORKING_DIR/${cname}/${cname}_fingerprint.txt ]]
then

  if [[ "$t" == "H3K27ac" || "$t" == "H3K4me3" || "$t" == "H2AFZ" || "$t" == "H3ac" || "$t" == "H3K4me2" || "$t" == "H3K9ac" ]]
  then
    bin_size=200
  else
    bin_size=1000
  fi

  echo "Experiment type: $t and bin size: $bin_size" >&2

  $DEEPTOOLS_2501/plotFingerprint -b $WORKING_DIR/${cname}/${cname}_dedup.bam $WORKING_DIR/${iname}/${iname}_dedup.bam -bs ${bin_size} -l ${cname} ${iname} --JSDsample $WORKING_DIR/${iname}/${iname}_dedup.bam --outQualityMetrics $WORKING_DIR/${cname}/${cname}_fingerprint.txt -plot $WORKING_DIR/${cname}/${cname}_fingerprint.png -p $n

  echo "return[run deeptools plotFingerprint]:$? "
else
  echo "deeptools plotFingerprint .txt exists: ${cname}_fingerprint.txt"
fi

if [[ "$t" == "Input" && "${cname}" == "${iname}" ]]
then
  js_dist=0
  chance_div=0
else
  js_dist=`grep ${cname} $WORKING_DIR/${cname}/${cname}_fingerprint.txt | cut -f 8`
  chance_div=`grep ${cname} $WORKING_DIR/${cname}/${cname}_fingerprint.txt | cut -f 12`
fi

#4.     Calculating FRiP scores
if [[ "$t" == "Input" ]]
then
  frip=0
else
  reads_under_peaks=`${SAMTOOLS_131}/samtools view -c -L ${p} $WORKING_DIR/${cname}/${cname}_dedup.bam`
  frip=$(echo "${reads_under_peaks}/${final_reads}" | bc -l)
fi

printf "ChIP_name\tInput_name\ttotal_reads\tmapped_reads\tdupped_reads\tdup_rate\tsingletons\tfinal_reads\tjs_dist\tchance_div\tfrip\n" > $WORKING_DIR/${cname}/${cname}_read_stats.txt
printf "%s\t%s\t%d\t%d\t%d\t%.4f\t%d\t%d\t%.4f\t%.4f\t%.4f\n" "$cname" "$iname" "$total_reads" "$mapped_reads" "$dupped_reads" "$dup_rate" "$singletons" "$final_reads" "$js_dist" "$chance_div" "$frip" >> $WORKING_DIR/${cname}/${cname}_read_stats.txt


