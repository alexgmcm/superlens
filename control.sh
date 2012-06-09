#!/bin/bash

#control.sh

#This script accepts the parameters for the metamaterials simulations as arguments
#and edits the fortran and matlab scripts to have those parameters

#It then compiles the fortran code, uses mv to rename the data files to a convenient name
#and then executes the matlab script to generate the plots

#This script will not loop over ranges of parameters so that it will be easy to call it to get one set
#but it is instead intended that this code can be executed by another loop script to generate ranges of parameters

#Expected usage:
#./control.sh eps2 mu2 eta thetai eps1 mu1 
#./control.sh 1.0 1.0 -2.0 -1.0 3 PI/4.0 
echo "Editing Fortran parameters to eps2= $3 , mu2= $4 , eta= $5 , thetai= $6 , eps1= $1 , mu1= $2"
eps2="$3"
mu2="$4"
eta="$5"
thetai="$6"
eps1="$1"
mu1="$2"

sed -e "s/^eps1=.*/eps1=$1/" congaussian1fixed.f90

echo "Completed."