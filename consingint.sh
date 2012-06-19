#!/bin/bash

#control.sh

#This script accepts the parameters for the metamaterials simulations as arguments
#and edits the fortran and matlab scripts to have those parameters

#It then compiles the fortran code, uses mv to rename the data files to a convenient name
#and then executes the matlab script to generate the plots

#This script will not loop over ranges of parameters so that it will be easy to call it to get one set
#but it is instead intended that this code can be executed by another loop script to generate ranges of parameters

#Expected usage:
#./control.sh eps2 mu2 eta thetai eps1 mu1 g
#./control.sh 1.0 1.0 5.0 1.0 3 "PI\/4.0" 3.0
echo "Editing Fortran parameters to eps1= $1 , mu1= $2 , eps2= $3 , mu2= $4 , eta= $5 , thetai= $6 , g=$7"
eps1="$1"
mu1="$2"
eps2="$3"
mu2="$4"
eta="$5"
thetai="$6"
g="$7"

#gets filename to put in for dat file saving/loading

striptheta=`echo $thetai | sed 's%\/%d%' `
#echo $striptheta
fname="'data/$striptheta rads $eta eta $g g $eps2 eps2gaussdielecfieldmap.dat'"
stripfname=`echo -n $fname | sed 's% %%g' ` 
#echo $stripfname

sed -e "s/^eps1=.*/eps1=$1/" \
-e "s/^mu1=.*/mu1=$2/" \
-e "s/^eps2=.*/eps2=($3,0.0)/" \
-e "s/^mu2=.*/mu2=$4/" \
-e "s/^eta=.*/eta=$5/" \
-e "s/^thetai=.*/thetai=$6/" \
-e "s/^ti='.*/ti='$6'/" \
-e "s/^g=.*/g=$7/" \
-e "s%write(filename,20).*%write(filename,20) $stripfname%" \
congaussian1fixed.f90 > congaussian1fixed.tmp
cp congaussian1fixed.tmp congaussian1fixed.f90
echo "Done."

echo "Editing MATLAB parameters..."

#get filename for plots
pname="$striptheta rads $eta eta $g g $eps2 eps2gaussdielecfieldmap.png');"
#echo $pname
strippname=`echo $pname | sed 's% %%g' `
#echo $strippname

sed -e "s/^eps2=.*/eps2=$3;/" \
-e "s/^mu2=.*/mu2=$4;/" \
-e "s/^eta=.*/eta=$5;/" \
-e "s/^thetai=.*/thetai='$6';/" \
-e "s/^g=.*/g=[$7];/" \
-e "s%load(.*%load($stripfname);%" \
-e "s%plots/im_.*%plots/im_$strippname%" \
-e "s%plots/real_.*%plots/real_$strippname%" \
congaussscript.m > congaussscript.tmp
cp congaussscript.tmp congaussscript.m 
echo "Done."

echo "Compiling Fortran code..."
gfortran congaussian1fixed.f90 -o a.out &>/dev/null
echo "Done."

echo "Executing Fortran code..."
./a.out &>/dev/null
echo "Done."

echo "Executing MATLAB code..."
matlab -nodisplay -nodesktop -r "congaussscript; quit" &>/dev/null
echo "Done."

#fix matlab messing with terminal
stty echo

echo "Completed."