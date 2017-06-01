#!/bin/bash

if [ $# -ne 2 ]; then
	echo "usage:">&2
	echo "extract.sh __file__ _line_">&2
	exit 1
fi


echo `head -n $2 $1 | tail -n 1`
