#!/bin/sh

seq 90 -10 10 > maxthetalist.txt
seq 3.0 1.0 3.0 > secondinterfacelist.txt
seq 30 2 30 > combcutofflist.txt


while read combcutoff
	do
	while read secint
		do
		echo "$secint"
		while read maxtheta           
			do   
			echo "$maxtheta"        
			./singprop.out $secint $maxtheta $combcutoff
			./matlab_batcher4.sh singpropscript $secint $maxtheta $combcutoff

		done <maxthetalist.txt
		#./matlab_batcher.sh resolutionplotscliscript $secint
	done <secondinterfacelist.txt
done <combcutofflist.txt




#rm list.txt