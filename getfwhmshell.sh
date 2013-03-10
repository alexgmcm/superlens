#!/bin/sh

seq 90 -10 10 > maxthetalist.txt
while read maxtheta           
			do   
			echo "$maxtheta"        
			./mastertest.out 3 $maxtheta 30 1 1 0.001
			./matlab_batcher4.sh fwhmioscript $maxtheta 0.001 \'combined\'

		done <maxthetalist.txt
