from utilsm import *
import os
import sys



def findfiles(base, pattern):
	cmd = "find {0} -name '{1}'".format(base, pattern)
	print cmd
	p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	return [e.strip() for e in p.stdout.readlines()]


def byino(fs):
	hashed = dict()
	for e in fs:
		ino = os.stat(e).st_ino 
		if not ino in hashed: hashed[ino] = list()
		hashed[ino].append(e)
	hashed2 = dict()
	for k,v in hashed.items():
		hashed2[k] =  sorted(v, key = lambda x: os.path.basename(x))[0]
	return hashed2

def md5script(hashed):
	def cmd(f):
		if f.strip().endswith('bam'):
			return 'echo "{1} $(./headlessbam_md5 {0})"'.format(f, os.path.basename(f))
		else:
			return 'echo "{1} $(md5sum {0})"'.format(f, os.path.basename(f))
			
	return [cmd(v) for v in sorted(hashed.values(), key= lambda x: os.path.basename(x))]

def trackoutput(base, i):
	logerr('# looking in {0}\n'.format(base))
	bams = findfiles(base, '*.bam')
	narrowpeaks = findfiles(base, '*narrow*gz')
	bamsbyino = byino(bams)
	peaksbyino = byino(narrowpeaks)
	print writef('./computemd5s_{0}'.format(i),  ['#!/bin/bash'] + md5script(bamsbyino) + md5script(peaksbyino)) 
	print findfiles(base, 'qc.html')	



def main(args):
	targets = ['.'] if len(args) == 0 else args
	for i, arg in enumerate(targets):
		trackoutput(arg, i)



if __name__ == '__main__':
	main(sys.argv[1:])



