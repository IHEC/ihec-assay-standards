from utilsm import * 
import sys
import glob

def check(data):
	ok = data['match_overall'] and len(data["failed_task_labels"]) == 0
	return ok


def main(args):
	N = 14
	home = args[0]
	n_ok = 0
	logs = glob.glob('{0}/*json'.format(home))
	for log in logs:
		data = jloadf(log)
		if not check(data):
			print '# failed: ' + log, data['match_overall'], data["failed_task_labels"]
		else:
			print '# ok:' + log
			n_ok +=1
	
	status = {
		"#ok": n_ok,
		"#expected": N,
		"#failed": N - n_ok
	}
	print jsonp(status)

if __name__ == '__main__':
	main(sys.argv[1:])
