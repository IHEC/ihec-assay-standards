# ENCODE ChIP-Seq pipeline singularity wrapper

## ENCODE ChIP-Seq pipeline

See https://github.com/ENCODE-DCC/chip-seq-pipeline2/


## Running tests

First run `./get_encode_resources.sh` to get encode test dataset and hg38 genome files. 

By default it will use git over http. If you want to use ssh, then pass `ssh` as first argumnet

Then run `python chip.py -pullimage -bindpwd` . Bind pwd will mount the current directory (equivalent to arguments `-B $PWD`). 

This will write:

* piperunner.sh

* testrun_tasks.sh

* singularity_encode_test_tasks.sh

* singularity_wrapper.sh

If you are running in `Local` mode using using `./chip.py -pullimage -bindpwd $PWD/data_b $PWD/data_a` will mount `$PWD/data_b` as `/mnt/ext_1`, `$PWD/data_a` as `/mnt/ext_2` and so on. It binds `$PWD` to `/mnt/ext_0`.  

If you are running in `sge_singularity` or `singularity`, binding local paths to container paths where they are not both the same will not work. If this use case is unavoidable, then refer to ENCODE documentation on how to use the pipeline. This IHEC wrapper does not support this.  

This will also create the singularity image in `./images` .

Do `chmod +x ./*sh`

You can pass `-nobuild` if you hust want to regenerate the wrapper scripts without pulling the singularity image again. 

Check singularity version with `singularity --version` to make sure it's at least `2.5.2` .

To run ENCODE test tasks, do `./singularity_test_tasks.sh Local` . The output of tests will be written in ` test_tasks_results` . Make sure all test pass, by looking through jsons generated. 

Run `chip.py -get` to get IHEC ChIP test data for MCF10A cell line. 

Doing `python chip.py  -maketests` will write ChIP test configurations:

* ./v2/ihec/cemt0007_h3k4me3.json

* ./v2/ihec/cemt0007_h3k27me3.json

IHEC tests can be run with:

`./singularity_wrapper.sh ./v2/ihec/cemt0007_h3k4me3.json` and `./singularity_wrapper.sh ./v2/ihec/cemt0007_h3k27me3.json` 

The provided configuration files are for 75bp PET only. Standard configration files for SET and read lengths will be provided. Currently the only local mode is supported for singularity. Slurm support is on the way. The ENCODE documentation addresses both. 


# To do

These items are pending for calling this done - 

* Waiting on encode v2 image

* Adding md5 checks for v2 output. 








