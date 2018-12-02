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

To run ENCODE test tasks, do `singularity_encode_test_tasks.sh Local try1` . The output of tests will be written in `test_tasks_results_try1` . Make sure all test pass, by looking through jsons generated. The first argument is the confi argument to cromwell (see ENCODE pipeline documentation). Only Local is currently supported.

Run `chip.py -get` to get IHEC ChIP test data for MCF10A cell line. 

Doing `python chip.py  -maketests` will write ChIP test configurations:

* ./v2/ihec/cemt0007_h3k4me3.json

* ./v2/ihec/cemt0007_h3k27me3.json

IHEC tests can be run with:

`./singularity_wrapper.sh ./v2/ihec/cemt0007_h3k4me3.json` and `./singularity_wrapper.sh ./v2/ihec/cemt0007_h3k27me3.json` 

The provided configuration files are for 75bp PET only. Standard configration files for SET and read lengths will be provided. Currently the only local mode is supported for singularity. Slurm support is on the way. The ENCODE documentation addresses both. 


To compute md5s of generated file, use `trackoutput.py <output_dir_1> ...`. This will locate peak calls and bam files, and generate scripts `computemd5s_$i` to compute the md5s. Note the bam md5s are generated without teh bam header as that may contain full paths names. 


As an example, supose output of `./singularity_wrapper.sh ./v2/ihec/cemt0007_h3k4me3.json` is in `$PWD/h3k4me3_out`. So do 

    python trackoutput.py $PWD/h3k4me3_out
	chmod +x ./computemd5s_0
	./computemd5s_0 > log_h3k4me3
	python status_cemt.py log_h3k4me3 expected_md5s_h3k4me3.json 

This will match md5s for cemt0007 H3K4me3 analysis. And similarly for H3K27me3. 

    $ python status_cemt.py computemd5s_0.out ./expected_md5s_h3k27me3.json 
    ok ChIP-Seq.IX1239-A26688-GGCTAC.134224.D2B0LACXX.2.1.merged.nodup.pr2_x_ctl_for_rep1.pval0.01.500K.narrowPeak.gz 1c9554fe8b67e61fd7c69a1881ec2e3a
    ok conservative_peak.narrowPeak.hammock.gz b78724bb667cc7bbfece8a587c10c915
    ok ChIP-Seq.IX1239-A26688-GGCTAC.134224.D2B0LACXX.2.1.merged.nodup.pr1_x_ctl_for_rep1.pval0.01.500K.bfilt.narrowPeak.hammock.gz defd886ab7923b952e04ee033a722fac
    ok optimal_peak.narrowPeak.hammock.gz b78724bb667cc7bbfece8a587c10c915
    ok rep1-pr.overlap.bfilt.narrowPeak.hammock.gz b78724bb667cc7bbfece8a587c10c915
    ok rep1-pr.overlap.narrowPeak.gz a896c1ec4693ddbd2e098ffa901c1f2a
    ok optimal_peak.narrowPeak.gz 49fdef6c06796ab06e8ac2a1b88075d1
    ok rep1-pr.overlap.bfilt.narrowPeak.gz 49fdef6c06796ab06e8ac2a1b88075d1
    ok ChIP-Seq.IX1239-A26688-GGCTAC.134224.D2B0LACXX.2.1.merged.nodup.pr1_x_ctl_for_rep1.pval0.01.500K.narrowPeak.gz b1ae4fb3f2b68b3c8346c57fa04f476f
    ok ChIP-Seq.IX1239-A26688-GGCTAC.134224.D2B0LACXX.2.1.merged.nodup_x_ctl_for_rep1.pval0.01.500K.bfilt.narrowPeak.gz 55de2037c6657d1027fb6b625822fa8b
    ok ChIP-Seq.IX1239-A26688-GGCTAC.134224.D2B0LACXX.2.1.merged.nodup.pr1_x_ctl_for_rep1.pval0.01.500K.bfilt.narrowPeak.gz 018ad8f5f3158534320ed359563878d3
    ok ChIP-Seq.IX1239-A26688-GGCTAC.134224.D2B0LACXX.2.1.merged.nodup_x_ctl_for_rep1.pval0.01.500K.bfilt.narrowPeak.hammock.gz bf5cd1f743325a0161d4ab77b00af829
    ok ChIP-Seq.IX1239-A26688-GGCTAC.134224.D2B0LACXX.2.1.merged.nodup_x_ctl_for_rep1.pval0.01.500K.narrowPeak.gz 7a52f55148b47e2a48fac330e3672c96
    ok conservative_peak.narrowPeak.gz 49fdef6c06796ab06e8ac2a1b88075d1
    ok ChIP-Seq.IX1239-A26688-GGCTAC.134224.D2B0LACXX.2.1.merged.nodup.pr2_x_ctl_for_rep1.pval0.01.500K.bfilt.narrowPeak.gz 0f38658b68706ec12b5faded1141750e
    ok ChIP-Seq.IX1239-A26688-GGCTAC.134224.D2B0LACXX.2.1.merged.nodup.pr2_x_ctl_for_rep1.pval0.01.500K.bfilt.narrowPeak.hammock.gz b1ac6ab70d053b546f186080639252ed
    {
        "failures": 0
    }





# To do

These items are pending for calling this done - 

* Waiting on encode v2 image

* Adding md5 checks for v2 output. 








