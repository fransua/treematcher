#!/usr/bin/env python 


from ete3 import Tree
from sys import stdin

for line in stdin:
	print Tree(line)
