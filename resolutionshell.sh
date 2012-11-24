#!/bin/sh

seq 60 5 90 > maxthetalist.txt
seq 4.5 0.5 10.0 > secondinterfacelist.txt


while read maxtheta
do
	echo "$maxtheta"
while read secint           
do   
	echo "$secint"        
   ./test.out $secint $maxtheta
   ./matlab_batcher.sh resolutionsuperlenscliscript $secint $maxtheta
               
done <secondinterfacelist.txt
done <maxthetalist.txt





#rm list.txt