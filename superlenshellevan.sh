#!/bin/sh
seq 4.0 0.5 4.0 > list.txt

while read line           
do   
	echo "$line"        
   ./test.out $line
   ./matlab_batcher.sh superlenscliscriptevan $line
               
done <list.txt





#rm list.txt