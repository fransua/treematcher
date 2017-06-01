#!/bin/bash

if [ $# != 2 ]; then
	echo "give a file and a new file name">&2
elif [ ! -f $1 ]; then
	echo "give a file">&2
	exit 1
fi

new_file=$2

if [ -f $new_file ]; then
	rm $new_file
fi

while read line; do
	echo $line | cut -f 5 -d " " >> $new_file
done < "$1"


echo "Completed!"
