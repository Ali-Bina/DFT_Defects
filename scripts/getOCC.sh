#!/bin/bash

#extracts occupations, energies and weights from QE SCF output
#formatted as columns in the order:Kpoints    Energies     occupations      weights
#parameters: pw.x scf/relaxation output file = $1

Name=$1

#extracts the relevant portion of the scf output containing energies, occupations vs k

begLine=$(( `grep -n 'End of self-consistent calculation' $Name | cut -d':' -f 1 | tail -n 1` + 1 ))
endLine=$(( `grep -n '!' $Name | cut -d':' -f 1 | tail -n 1` - 3 ))
head -n $endLine $Name | tail -n $(( $endLine - $begLine )) > tmp.txt



flag=0


#Stores energies for a given k
declare -a Karr

#stores occupations for a given k
declare -a Oarr

 
#read lines one by one
while read -r line; do


    #once the nbnd energies and occupations corresponding to a single k have both been scraped
    #write them as columns to a file
    if (( ${#Karr[@]} > 0 && ${#Oarr[@]} > 0 && ${#Karr[@]} == ${#Oarr[@]}));then
	
	nband=${#Karr[@]}
	
	#appends energies and occupations as columns to file
	for (( i=0; i < ${#Karr[@]}; i++ ))
	    {
		
		echo -e $kvect"\t${Karr[$i]}\t${Oarr[$i]}" >> 'tmp2.txt'
		
	    }

	    #delete arrays
	    unset Karr
	    unset Oarr
	fi
    
    if [[ "$line" =~ 'k =' ]];then
     
	flag=1
	kvect="`echo $line | grep -o '[0-9][0-9]*\.[0-9]*'`"


    elif [[ "$line" =~ 'occupation numbers' ]];then
	
	flag=0
    fi
   
    #stores energies 
    if [[ $flag == 1 && "$line" =~ ^-*[0-9].* ]];then
	Karr+=($line)

	
    #stores occupations
    elif [[ $flag == 0 && "$line" =~ ^-*[0-9].* ]];then
	Oarr+=($line)
    fi
    
	
	done < "tmp.txt"

	    
rm tmp.txt

#extracting the weights


firstLine=`cat $1 | grep -n 'number of k points='`
numKs=`echo $firstLine | cut -d' ' -f 6`
begLine=`echo $firstLine | cut -d':' -f 1`


array=( `head -n $(( $begLine + $numKs + 1 )) $1 | tail -n $numKs | cut -d'=' -f 3 | tr -s [:blank:]` )

for elem in ${array[@]};do
    for (( j=0; j < ${nband}; j++ )){
	echo $elem >> tmp3.txt
    }
done


paste tmp2.txt tmp3.txt > tmp4.txt

rm tmp2.txt tmp3.txt

cat tmp4.txt | awk 'BEGIN {print  "Kvector\t\t\tenergy\toccupations\tweights"} {print $0}' 

rm tmp4.txt
