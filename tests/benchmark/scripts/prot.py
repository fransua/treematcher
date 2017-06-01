#!/usr/bin/env python

from ete3 import Tree
from treematcher import TreePattern
import cProfile

total = 0
correct = 0
error = 0 
matches = 0

profiler = cProfile.Profile()

file_name ="art_trees.txt"

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
			#print tree
			#print pt
			res = pt.find_match(tree, maxhits=None)
			profiler.disable()
			correct += 1
			if len(list(res)) > 0:
				matches += 1
		except:
			error += 1
		#print str(total) 
			
						
		#except:
			#traceback.print_exception() > errors_log
			#error += 1
	print file
	print "total: " + str(total)
	print "succes: " + str(correct)
	print "error : " + str(error)
	print "rate  : " + str ( float(correct)/float(total)*100)
	print "match : " + str(matches) + " / " + str(correct)	

	profiler.print_stats(sort="time")





