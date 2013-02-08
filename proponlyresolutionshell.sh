#!/bin/sh

seq 90 -10 90 > maxthetalist.txt
seq 3.0 1.0 3.0 > secondinterfacelist.txt


while read secint
	do
	echo "$secint"
	while read maxtheta           
		do   
		echo "$maxtheta"        
		./proptest.out $secint $maxtheta
		./matlab_batcher3.sh proponlyresolutionsuperlenscliscript $secint $maxtheta

	done <maxthetalist.txt
	#./matlab_batcher.sh resolutionplotscliscript $secint
done <secondinterfacelist.txt





#rm list.txt