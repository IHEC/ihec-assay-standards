#!/bin/bash


BASE=$1
BACKEND=$2
H=$BASE
tag=${3:-""}

chmod +x $H/testrun_tasks.sh
testsOut=$H/test_tasks_results_"$tag"
mkdir $testsOut || true
cd $BASE/chip-seq-pipeline2/test/test_task
echo "__container__:$BASE,$BACKEND,$PWD"

for t in test_bam2ta test_bwa  test_choose_ctl test_filter test_idr test_macs2 test_merge_fastq test_overlap test_pool_ta test_reproducibility test_spp test_spr test_trim_fastq test_xcor; do
#for t in test_bam2ta; do
  echo "# started: $t $(date)"
  $H/testrun_tasks.sh $PWD/$t.wdl $PWD/$t.json $testsOut/$t.test_task_output.json $BACKEND
  echo "# end: $t $(date) $?"
  echo "ok___________________"
done



