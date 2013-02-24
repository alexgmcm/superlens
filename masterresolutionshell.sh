#!/bin/sh
#need to write matlabscript
seq 90 -10 90 > maxthetalist.txt
seq 3.0 1.0 3.0 > secondinterfacelist.txt
seq 30 2 30 > combcutofflist.txt
MODEFLAG=0
IMAGEFLAG=0

while read combcutoff
	do
	while read secint
		do
		echo "$secint"
		while read maxtheta           
			do   
			echo "$maxtheta"        
			./mastertest.out $secint $maxtheta $combcutoff $MODEFLAG $IMAGEFLAG
			./matlab_batcher4.sh combinedresolutionsuperlenscliscript $secint $maxtheta $combcutoff

		done <maxthetalist.txt
		#./matlab_batcher.sh resolutionplotscliscript $secint
	done <secondinterfacelist.txt
done <combcutofflist.txt




#rm list.txt