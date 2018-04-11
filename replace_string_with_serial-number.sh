#!/bin/bash
######################################################
## determine the number of occurences of $string
## replace each single occurence with a serial number
## starting from one until #$string is reached
######################################################

target=assembly.fa
string="gi.*"			## string starting with "gi" followed by arbitrary characters/numbers
j=`grep "$string" assembly.fa | wc | awk '{print $1}'`

for ((i=1; i<=$j; i++))
do
	## sed -inline "start-from-line-0/find-$string/replace/$string/with-serial-number/" target-file
	sed -i "0,/$string/s/$string/$i/" assembly.fa
	echo "... replaced $i of $j"
done
