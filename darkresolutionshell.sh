#!/bin/sh

seq 90 -10 90 > maxthetalist.txt
seq 3.0 1.0 3.0 > secondinterfacelist.txt
seq 3 1 10 > darkcutofflist.txt


while read darkcutoff
do
	echo "$darkcutoff"
	while read secint
		do
		echo "$secint"
		while read maxtheta           
			do   
			echo "$maxtheta"        
			./darktest.out $secint $maxtheta $darkcutoff
			./matlab_batcher4.sh darkresolutionsuperlenscliscript $secint $maxtheta $darkcutoff

		done <maxthetalist.txt
		#./matlab_batcher.sh resolutionplotscliscript $secint
	done <secondinterfacelist.txt
done <darkcutofflist.txt




#rm list.txt