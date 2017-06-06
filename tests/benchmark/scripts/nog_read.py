#!/usr/bin/env python

import sys
import time
from ete3 import Tree


if len(sys.argv) < 1:
	print "File parameter not found"
	sys.exit(1)




lines = 0
errors = 0
total = 0
start = time.time()
with open(sys.argv[1]) as file:
	for line in file:
		lines += 1
		if len(line) > 1:
			try:
				Tree(line)
				total += 1
			except:
				errors += 1
				print "error at " + str(lines)

end = time.time()	
completed_at = end - start 
print str(lines) + " lines found"
print "Success: " + str(total)
print "Erros  : " + str(errors)
print "Empty  : " + str(lines - (total + errors))
print "Script completed in " + str(completed_at/60) + " minutes"
