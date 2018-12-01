import json
import subprocess
import sys
import os
import glob
import shutil

def jsonp(m):
	return json.dumps(m, sort_keys=True, indent=4)


def jdumpf(f, m):
	with open(f, 'w') as outfile:
		outfile.write(jsonp(m))
	return f

def dumpf(f, data):
	with open(f, 'w') as outfile:
		outfile.write(data)
	return f

def writef(f, entries, sep = '\n'):
	with open(f, 'w') as outfile:
		for e in entries:
			outfile.write(e)
			outfile.write(sep)	
	return f
	
def linesf(f):
	with open(f) as infile:
		return infile.readlines()

	
def jloadf(f):
	with open(f) as infile:
		return json.load(infile)


def logerr(m):
	sys.stderr.write('{0}'.format(m))

def logerrn(m):
	sys.stderr.write('{0}\n'.format(m))


def shell(cmd, assert_ok=False):
	p = subprocess.Popen(cmd, shell=True)
	exit = p.wait()
	logerrn('#exit:{1} {0}'.format(cmd,exit))
	if assert_ok:
		assert exit == 0
	return exit == 0

def mkdirs(path):
	try:
		os.makedirs(path)
	except OSError as err:
		logerr(str(err) + '\n')
	
def movefile(s, d):
	shutil.move(s, d)




def by_keyvalue(alist, k, v):
	hashed = dict()
	for e in alist:
		ke = k(e)
		ve = v(e)
		if not ke in hashed:
			hashed[ke] = list()
		hashed[ke].append(ve)
	return hashed




