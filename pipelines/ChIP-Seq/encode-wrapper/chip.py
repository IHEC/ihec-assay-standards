import json
import subprocess
import sys
import os
import glob
import shutil



from utilsm import *



debug_mode = False

def base():
	return os.path.dirname(os.path.realpath(__file__))


def wget(url, debug=debug_mode):
	logerr('getting: {}\n'.format(url))
	if debug:
		logerr(' ..debug: wget {0}\n'.format(url))
		dumpf(os.path.basename(url), 'test:{0}'.format(url))
		return
	p = subprocess.Popen('wget ' + url ,shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	#p = subprocess.Popen(['wget', url, '--directory-prefix', './test_data'] ,shell=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	for line in p.stdout.readlines():
		logerr(line)
	return p.wait()
	
	
	
def get_hg38_resources(home):
	base = os.path.abspath(os.getcwd())
	mkdirs('hg38_resources/genome_hg38/bwa_index')
	os.chdir('./hg38_resources')
	for f in [
			'http://www.epigenomes.ca/data/CEMT/resources/chip_v2/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta',
			'http://www.epigenomes.ca/data/CEMT/resources/chip_v2/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.tar',
			'http://www.epigenomes.ca/data/CEMT/resources/chip_v2/hg38.blacklist.bed.gz',
			'http://www.epigenomes.ca/data/CEMT/resources/chip_v2/hg38.chrom.sizes',
		]:
		wget(f)
	movefile('hg38.blacklist.bed.gz', './genome_hg38/')
	movefile('GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta', './genome_hg38/')
	movefile('hg38.chrom.sizes', './genome_hg38/')
	movefile('GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.tar', './genome_hg38/bwa_index/')
	config = '\n'.join([
		"blacklist {0}/genome_hg38/hg38.blacklist.bed.gz",
		"chrsz     {0}/genome_hg38/hg38.chrom.sizes",
		"gensz     hs",
		"bowtie2_idx_tar   /dev/null",
		"bwa_idx_tar       {0}/genome_hg38/bwa_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.tar",
		"ref_fa    {0}/genome_hg38/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta",
	]).format(base) + '\n'
	logerr( dumpf('./hg38_local.tsv', config) + '\n' )
	base_config = jdumpf('./base_config.json', { 'chip.genome_tsv' : os.path.abspath('./hg38_local.tsv'), 'base' : base  })
	os.chdir(home)
	return base_config


def existing_ref_config(configfile):
	home = os.path.abspath(os.getcwd())
	config = jloadf(configfile)
	hashed = { k : os.path.realpath(config[k]) if k in ['blacklist', 'chrsz', 'ref_fa', 'bwa_idx_tar'] else config[k] for k in config }
	mkdirs('./hg38_resources')
	os.chdir('./hg38_resources')
	config = '\n'.join([e.strip() for e in '''blacklist	{blacklist}
				chrsz	{chrsz}
				gensz	{gensz}
				bowtie2_idx_tar	{bowtie2_idx_tar}
				bwa_idx_tar	{bwa_idx_tar}
				ref_fa	{ref_fa}'''.format(**hashed).splitlines()  ]) + '\n'
	logerr( './hg38_resources/' + dumpf('./hg38_local.tsv', config) + '\n' )
	base_config = jdumpf('./base_config.json', { 'chip.genome_tsv' : os.path.abspath('./hg38_local.tsv'), 'base' : home  })
	os.chdir(home)
	return home  + '/hg38_resources/base_config.json'


def rm(target):
	try:
		shutil.rmtree(target)
	except OSError as e:
		logerr("# error: {0} / {1}".format(target, e.strerror))


def get_test_data(configfile, home):
	config = jloadf(configfile)
	os.chdir('./v2/ihec/test_data')
	oks = dict()
	for k in config['data']:
		oks[k] = False
		for url in config['data'][k]:
			if wget(url) == 0:
				oks[k] = True
				break
			else:
				logerr('# failed downloading:' + url)
				incomplete = glob.glob('./' + os.path.basename(url))
				if len(incomplete) > 0:
					assert len(incomplete) == 1, incomplete
					shutil.remove(incomplete[0])
					logerr('# removed failed download.. ' + incomplete[0])
	os.chdir(home)	
	for k in oks:
		assert oks[k], ['could not download all test data', k]



def make_tests():
	mcf10a = ['cemt0007_h3k4me3_template.json', 'cemt0007_h3k27me3_template.json'] 
	
	if os.path.isfile('./hg38_resources/base_config.json'):
		config = jloadf('./hg38_resources/base_config.json')
		base = config['base']
	else:
		base = os.path.abspath(os.getcwd())

	def fix(fname, base):
		assert fname.endswith('_template.json')
		out = './v2/ihec/{0}.json'.format(fname[0:-len('_template.json')])
		config = jloadf(fname)
		return dumpf(out, jsonp(config).replace('{0}', base))
	
	for f in mcf10a:
		print 'written:', fix(f, base)

def write_testrun(config):
	with open('testrun_template.sh') as infile:
		logerr(dumpf('{0}/testrun.sh'.format(config['home']),  infile.read().format(**config)) + '\n')
	with open('testrun_tasks_template.sh') as infile:
		logerr(dumpf('{0}/testrun_tasks.sh'.format(config['home']),  infile.read().format(**config)) + '\n')
	
	logerrn(dumpf('./singularity_test_tasks.sh', '#!/bin/bash\n\necho "home:$PWD"\nwhich singularity\nsingularity exec {container_image} {home}/encode_test_tasks_run.sh\n\n'.format(**config)))
	return dumpf('./singularity_test.sh', '#!/bin/bash\n\necho "home:$PWD"\nwhich singularity\nsingularity exec {container_image} {home}/testrun.sh $1\n\n'.format(**config))



def singularity_pull_image(home, config, debug=debug_mode):
	imageurl = 'docker://quay.io/encode-dcc/chip-seq-pipeline:v2'
	imageurl = 'docker://quay.io/encode-dcc/chip-seq-pipeline:v1.1'
	os.chdir('./images')
	if debug:
		dumpf('./debug.img', 'test:{0}'.format('singularity'))
	else:
		cmd = 'singularity pull {0}'.format(imageurl)
		logerr('# .. ' +  cmd + '\n')
		if not '-nobuild' in config:
			shell(cmd, assert_ok = True)
	
	images = glob.glob('./*img')
	assert len(images) == 1
	image_label = 'chip_seq_pipeline_v_1_1'
	image_name = '{0}.simg'.format(image_label)
	logerr('# pulled image: {0}, moved: {1}\n'.format(images[0], image_name))
	os.rename(images[0], image_name)
	image_path = os.path.abspath(image_name)
	os.chdir(home)
	container = jdumpf('./v2/singularity_container.json',  {
		"default_runtime_attributes" : {
			"singularity_container" : image_path,
			"singularity_instance_name": image_label
		}
	})
		
	shell('singularity exec {0} cp /software/chip-seq-pipeline/chip.wdl ./v2'.format(image_path), assert_ok=True)
	logerr('# copied /software/chip-seq-pipeline/chip.wdl to ./v2/chip.wdl\n')
	return {
		"container_image":image_path,
		"home" : home,
		"container" :  os.path.abspath(container),
		"wdl" : "{0}/v2/chip.wdl".format(home),
		"backend" : "{0}/backend.conf".format(home)
	}
def main(args):
	home = base()
	logerr('# prefix {0}\n'.format(home))
	mkdirs('./hg38_resources')
	mkdirs('./images')
	mkdirs('./v2/ihec/test_data')
	
	if '-clean' in args:
		for d in ['./v2', './images', './hg38_resources']:
			logerr('# removing {0}\n'.format(d))
			rm(d)
		logerr('rm -rf ./v2/ images/ hg38_resources/ \n')
	
	if '-getref' in args:
		get_hg38_resources(home)
	
	if '-refconfig' in args:
		logerr(existing_ref_config('./ref_config.json') + '\n')	

	if '-get' in args:
		get_test_data('./test_config.json', home)
	
	if '-pullimage' in args:
		container_config = singularity_pull_image(home, args, debug = False)
		container = write_testrun(container_config)
		logerr('# container: {0}\n'.format(container))
	
	if '-maketests' in args:
		make_tests()
	
	
	
	
	
if __name__ == '__main__':
	main(sys.argv[1:])
