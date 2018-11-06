# ENCODE ChIP-Seq pipeline singularity wrapper

## ENCODE ChIP-Seq pipeline

See https://github.com/ENCODE-DCC/chip-seq-pipeline2/


## Running tests

First run `./get_encode_resources.sh` to get encode test dataset and hg38 genome files. 

By default it will use git over http. If you want to use ssh, then pass `ssh` as first argumnet

Then run `python chip.py -pullimage` . This will write:

* testrun.sh

* testrun_tasks.sh

* singularity_test_tasks.sh

* singularity_test.sh


As well as creating the singularity image in `./images` .

Do `chmod +x ./*sh`

Check singularity version with `singularity --version` to make sure it's at least `2.5.2` .

To run ENCODE test tasks, do `./singularity_test_tasks.sh` . The output of tests will be written in ` test_tasks_results` . Make sure all test pass, by looking through jsons generated. 

Run `chip.py -get` to get IHEC ChIP test data for MCF10A cell line. 

Doing `python chip.py  -maketests` will write ChIP test configurations:

* ./v2/ihec/cemt0007_h3k4me3.json

* ./v2/ihec/cemt0007_h3k27me3.json

IHEC tests can be run with:

`./singularity_test.sh ./v2/ihec/cemt0007_h3k4me3.json` and `./singularity_test.sh ./v2/ihec/cemt0007_h3k27me3.json` 


The provided configuration files are for 75bp PET only. Standard configration files for SET and read lengths will be provided. Currently the only local mode is supported for singularity. Slurm support is on the way. The ENCODE documentation addresses both. 


# To do

These items are pending for calling this done - 

* Waiting on encode v2 image

* Adding md5 checks for v2 output. 








