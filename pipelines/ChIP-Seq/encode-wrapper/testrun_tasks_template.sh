#!/bin/bash

unset PYTHONPATH
unset R_LIBS_USER
which R
which python
which java

echo "paths: $R_LIBS_USER $PYTHONPATH"
echo $PATH


CROMWELL_HOME="{home_mnt}"
BACKEND_CONF="{backend}"
WORKFLOW_OPT="{container}"
BACKEND="$4"

WDL="$1"

PREFIX=$(basename $WDL .wdl)
METADATA="$PREFIX".metadata.json # metadata
RESULT=$3


jobFile=$2 

java -jar -Dconfig.file=$BACKEND_CONF -Dbackend.default=$BACKEND $CROMWELL_HOME/cromwell-34.jar run $WDL -i $jobFile -o $WORKFLOW_OPT -m $METADATA
echo "return:$?"

cat $METADATA | python -c "import json,sys;obj=json.load(sys.stdin);print(obj['outputs']['"$PREFIX".compare_md5sum.json_str'])" > $RESULT
cat $RESULT





