#!/bin/bash

#user_input() {
#	read answer
#
#	if [ ! answer]; then answer=$deafult_ans; fi
#
#	if [ $answer == "y"]; then	pass=1
#	elif [ $answer == "Y"]; then pass=1
#	else	pass=0
#	fi
#}


#if [ $# -eq 0]; then
#	if [ ! -d ./tress ]; then
#		echo "trees folder does not exist or it is empty do you want to download"
#		echo -n "a default set of trees? (y/n | y): "
#		user_input
#
#
#
#
#
#if [ $# -eq 0 ]; then
#		echo "No file is passed as parameter."
#		echo "All files will be scaned. It may take a lot of time."
#		echo -n "Do you want to continue? (y/n) "
#
#
#fi

if [ ! -d prfl ]; then
	mkdir prfl
fi

count=0
total=$#

for cFile in "$@"; do
	cFile=`echo $cFile | sed -r "s/\//\\\\\\\\\//g"`
	sed -r "s/__file__name__/$cFile/" ./scripts/prototype.py > ./scripts/temp_py.py
	name=`basename $cFile`
	python ./scripts/temp_py.py > prfl/prfl_${name}
	rm ./scripts/temp_py.py
	count=`expr $count + 1`
	echo "${count}/${total}: $name - completed!"

done

./results.sh
