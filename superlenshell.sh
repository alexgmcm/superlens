#!/bin/sh
seq 0.01 0.5 4.1 > list.txt

while read line           
do   
	echo "$line"        
   ./test.out $line
   ./matlab_batcher.sh superlenscliscript $line
               
done <list.txt





#rm list.txt