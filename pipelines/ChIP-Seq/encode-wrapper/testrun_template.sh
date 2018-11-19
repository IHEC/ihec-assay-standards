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
BACKEND=$2 #"{backend_default}"
CHIP="{wdl}"

jobFile=$1 

java -jar -Dconfig.file=$BACKEND_CONF -Dbackend.default=$BACKEND cromwell-34.jar run $CHIP -i $jobFile -o $WORKFLOW_OPT
echo "return:$?"
