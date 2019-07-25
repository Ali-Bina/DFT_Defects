#!/bin/bash

#Looks into the output directory from the current working directory
#Expects as input the prefix of files to be processed ie looks for files with names of the form $1_num.out
#Where num is the value of the relavent parameter (k grid, ecut, alat,...)


for output in ./outputs/${1}*.out
do

echo "`echo $output | sed "s/^\.[/]outputs[/]${1}_\([0-9]*\.*[0-9]*\)\.out$/\1/"` `grep ! $output | tr -s [:blank:] | cut -d' ' -f5`"
    
done
