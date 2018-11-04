from utilsm import *
import os
import sys



def findfiles(base, pattern):
	cmd = "find {0} -name '{1}'".format(base, pattern)
	print cmd
	p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	return [e.strip() for e in p.stdout.readlines()]


def byino(fs):
	hashed = {os.stat(e).st_ino : e  for e in fs}
	for k,v in hashed.items():
		print '#', os.path.basename(v)
	return hashed

def md5script(hashed):
	return ['echo "{1} $(md5sum {0})"'.format(v, os.path.basename(v)) for v in sorted(hashed.values(), key= lambda x: os.path.basename(x))]

def main(args):
	base = '.' if len(args) == 0 else args[0]
	logerr('# looking in {0}\n'.format(base))
	bams = findfiles(base, '*.bam')
	narrowpeaks = findfiles(base, '*narrow*gz')
	bamsbyino = byino(bams)
	peaksbyino = byino(narrowpeaks)
	print writef('./computemd5_2',  ['#!/bin/bash'] + md5script(bamsbyino) + md5script(peaksbyino)) 
	




if __name__ == '__main__':
	main(sys.argv[1:])



