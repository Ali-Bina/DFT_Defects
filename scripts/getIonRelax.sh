#!/bin/bash


#extract final atomic coordinate from a calculation='relax', provided the output file as the first argument, write to standard out
#if time is provided as the second argument then it will resubmit the job with the new coordinates for the specified time

#Program expects an input and an output directory ,storing the QE inputs/outputs, as child directories of the current working directory
#Stupid idiosyncracy: '.' in file name only before file extention allowed, '_' used to separated prefix from job type, i.e. names like Ge_Te0 can cause problems

output=''
outputFile=$1
readFile=$outputFile
prefix=`echo "$outputFile" | rev | cut -d'/' -f 1 | rev | cut -d'.' -f 1 | cut -d'_' -f 2`
inputFile=`echo "$outputFile" | rev | cut -d'/' -f 1 | rev | cut -d'.' -f 1 | sed 's/$/.in/' | sed 's/^/inputs\//'`
time=$2

if [[ ! -e $inputFile ]]; then

    echo "Error input file not found"
    exit
fi

start=`grep -n 'CELL_PARAMETERS'  "$inputFile" | cut -d':' -f 1`
cellParam=`head -n $(( $start + 3 ))  $inputFile  | tail -n 3`
output=`printf "$output\n$cellParam\n"`

#finds the first line of the final atomic coordinates
begin=`grep -n --text 'ATOMIC_POSITIONS' $outputFile | tail -n 1 | cut -d':' -f 1`

#if no new coordinates written to output then output input coordinates
if [[ $begin == '' ]]; then
    
    begin=`grep -n --text 'ATOMIC_POSITIONS' $inputFile | tail -n 1 | cut -d':' -f 1`
    readFile=$inputFile
fi

#Reads file upto begin the appends coordinates to file named $name
i=0
nAtoms=0
newCoords=''
while IFS= read -r line; do
    
    if (( $i >= $begin )) && [[ $line =~ ^[a-zA-Z][a-zA-Z].*[0-9]  ]]; then

	let nAtoms++

	if (( $nAtoms == 1  )); then
	    
	    newCoords=$line

	else

	    newCoords=`printf "$newCoords\n$line"`

	fi
	
    fi
    
    let i++
    
done < "$readFile"
output=`printf "$output\n$newCoords"`

#Metadata:
#final number of scf steps before crash
finalSCF=`grep 'number of scf cycles' $outputFile | tail -n 1`
#final number of bfgs steps before crash
finalBFGS=`grep 'number of bfgs steps' $outputFile | tail -n 1`
#Total CPU time before crash
cpuTime=`grep 'total cpu time spent up to now is' $outputFile | tail -n 1`

echo -e "$finalSCF\n$finalBFGS\n$cpuTime" > "${prefix}_metadata"

#resubmit job using new coordinates and information from previous input
if [[ "$time" == "" ]];then

    echo "$output"

else
    
    #copy over all Parameters from input up to ATOMIC_POSITIONS 
    i=0
    stop=`grep -n 'ATOMIC_POSITIONS' $inputFile | cut -d':' -f 1`
    while IFS= read -r line; do
	let i++
	submit=`echo -e  "$submit\n$line"`
	
	if (( $i >= $stop )); then
       
	    break
	
	fi 
	
	
    
    done < "$inputFile"

    #Add new Atomic positions from previous relaxation
    submit=`echo -e "$submit\n$newCoords"`

    #add the rest of input:
    cont=$(( $stop + $nAtoms ))
    fileLen=`cat $inputFile | wc -l`
    remainder=`tail -n $(( $fileLen - $cont )) $inputFile`

    #add the remaninder of the old input file to the new input file
    submit=`echo -e "$submit\n$remainder"`
    echo -e "$submit"
    echo -e "$submit" > "${inputFile}"

    echo "Submition Time: $time"
    echo "Job Prefix: $prefix"
    sbatch --time="$time" --job-name=${prefix} ./jobs/pwSubmit.sh "${inputFile}"
fi

