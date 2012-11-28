#!/bin/sh

seq 90 -10 10 > maxthetalist.txt
seq 3.0 1.0 3.0 > secondinterfacelist.txt


while read secint
	do
	echo "$secint"
	while read maxtheta           
		do   
		echo "$maxtheta"        
		./test.out $secint $maxtheta
		./matlab_batcher3.sh resolutionsuperlenscliscript $secint $maxtheta

	done <maxthetalist.txt
	./matlab_batcher.sh resolutionplotscliscript $secint
done <secondinterfacelist.txt





#rm list.txt