#!/bin/sh
#seq 1.5 1.2 > list.txt

while read line           
do   
	echo "$line"        
   ./test.out $line
   ./matlab_batcher.sh darksingintscript $line
               
done <darklist.txt





#rm list.txt