#!/usr/bin/env python

from ete3 import Tree
from treematcher import TreePattern
import cProfile
from sys import stderr

total = 0
correct = 0
error = 0
matches = 0

profiler = cProfile.Profile()

file_name ="__file__name__"

with open(file_name) as file:
	for line in file:
		try:
			total += 1
			tree = Tree(line)
			expanded = str(tree).splitlines()
			first = expanded[1].strip().replace('\\', '').replace('/', '')
			if first[0] == "-":
				first = first[1:]
			last  = expanded[len(expanded)-1].strip().replace('\\','')
			if last[0] == "-":
				last = last[1:]
			pt = "((" + str(first) + ")*), ((" + str(last) + ")*)@;"
			#pt = "((" + str(first) + ")*)@;"
			#pt = "((" + str(last) + ")*)@;"
			pt = TreePattern(pt)
			profiler.enable()
			res = pt.find_match(tree, maxhits=None)
			profiler.disable()
			correct += 1
			if len(list(res)) > 0:
				matches += 1
			else:
				print str(total)
		except:
			error += 1

		#except:
			#traceback.print_exception() > errors_log
			#error += 1
	print file.name
	print "total  : " + str(total)
	print "success: " + str(correct)
	print "error  : " + str(error)
	print "rate   : " + str ( float(correct)/float(total)*100) + "%"
	print "matches: " + str(matches) + " / " + str(correct)
	print "m_rate : " + str(float(matches) / float(correct) * 100) + "%"

	profiler.print_stats(sort="time")
