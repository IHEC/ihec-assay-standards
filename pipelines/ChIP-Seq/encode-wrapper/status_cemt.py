from utilsm import *






def main(args):
	def k(x):
		try:
			return  x.split()[0]
		except:
			return '-'
	def v(x):
		try:
			return x.split()[1]
		except:
			return '-'


	log = by_keyvalue(linesf(args[0]), k, v)
	expected = jloadf(args[1]) 
	n_fails = 0
	for k in expected:
		if not k in log:
			n_fails += 1
			print 'missing', k
		else:
			found = log[k]
			assert len(found) == 1, 'non unique md5: '+ str(k, found)
			if not  found[0] == expected[k][0]:
				n_fails += 1
			else:
				print 'ok', k, expected[k][0]
	
	result = {
		"failures":n_fails
	}
	print jsonp(result)


if __name__ == '__main__':	
	import sys
	main(sys.argv[1:])

