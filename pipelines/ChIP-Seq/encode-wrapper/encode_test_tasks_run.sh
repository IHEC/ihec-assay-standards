#!/bin/bash


BASE=$1
H=$BASE
chmod +x $H/testrun_tasks.sh
mkdir $H/test_tasks_results || true
cd $BASE/chip-seq-pipeline2/test/test_task
echo $BASE
for t in test_bam2ta test_bwa  test_choose_ctl test_filter test_idr test_macs2 test_merge_fastq test_overlap test_pool_ta test_reproducibility test_spp test_spr test_trim_fastq test_xcor; do
#for t in test_bam2ta; do
  echo "# started: $t $(date)"
  $H/testrun_tasks.sh $BASE/$t.wdl $BASE/$t.json $H/test_tasks_results/$t.test_task_output.json 
  echo "# end: $t $(date) $?"
  echo "ok___________________"
done



