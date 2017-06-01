#!/bin/bash

#echo -e "file\t\tlines\ttrees\tnot_trees\ttree_rate\tmatches\tmatches_rate"
printf "%-20s%-10s%-10s%-10s%-10s%-10s%-10s" "file" "lines" "trees" "not_trees" "tree_rate" "matches" "match_rate"
printf "\n"
for file in ./prfl/prfl_* ; do
	#name=`echo "$file" | cut -d"/" -f 2`
	name=`basename $file`

	total=`head -n 10 $file  | egrep "^total" | cut -d ":" -f 2`
  ntrees=`head -n 10 $file | egrep "^success"| cut -d ":" -f 2`
	errors=`head -n 10 $file | egrep "^error" | cut -d ":" -f 2`
	s_rate=`head -n 10 $file | egrep "^rate"  | cut -d ":" -f 2`
	matches=`head -n 10 $file| egrep "^matches" | cut -d ":" -f 2 | cut -d "/" -f 1`
	m_rate=`head -n 10 $file | egrep "^m_rate"| cut -d ":" -f 2`

	printf "%-20s%-10s%-10s%-10s%-10s%-10s%-10s\n" $name $total $ntrees $errors $s_rate $matches $m_rate

done


echo ""
