#!/bin/sh

matlab_exec=matlab
X="${1}(${2},${3})"
echo ${X} > matlab_command_${2}_${3}.m
cat matlab_command_${2}_${3}.m
${matlab_exec} -nojvm -nodisplay -nosplash < matlab_command_${2}_${3}.m
rm matlab_command_${2}_${3}.m


#Call it by entering:
#./matlab_batcher.sh myfunction myinput
#