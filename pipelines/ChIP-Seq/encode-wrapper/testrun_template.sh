#!/bin/bash

unset PYTHONPATH
unset R_LIBS_USER
which R
which python
which java

echo "paths: $R_LIBS_USER $PYTHONPATH"
echo $PATH

BACKEND_CONF="{backend}"
WORKFLOW_OPT="{container}"
BACKEND="Local"
CHIP="{wdl}"

jobFile=$1 #./v2/ihec_tests/cemt0007_h3k4me3.json

java -jar -Dconfig.file=$BACKEND_CONF -Dbackend.default=$BACKEND cromwell-34.jar run $CHIP -i $jobFile -o $WORKFLOW_OPT
echo "return:$?"
