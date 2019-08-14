#!/bin/bash

#generate a uniform grid of kpoints in coordinates of the reciprocal lattice
kpoints(){
for (( i=0; i < $1; i++)){
	for (( j=0; j < $1; j++)){
		for (( k=0; k < $1; k++)){
		        
		        awk -v i=$i -v j=$j -v k=$k -v n=$1 'BEGIN {printf "%0.15f\t%0.15f\t%0.15f\n", i/n, j/n, k/n}'
			
		    }
	    }
    }
    
}


#modifies Parameters of interest from a default job script
changeParams(){
    file=$1
    content=`cat $file`
    for param in ${modify[@]}; do

	val=`eval echo \\${$param}`
	
	if [[ $param == 'charge'  ]]; then
	    
	    content=`echo "$content" | sed -r  "s/tot_charge=[-]*[0-9]+/tot_charge=${val}/"`
	    
	elif [[ $param == 'ecut' ]];then

	    content=`echo "$content" | sed -r  "s/ecutwfc=[0-9]+/ecutwfc=${val}/"`

	elif [[ $param == 'lat' ]]; then

	    content=`echo "$content" | sed -r  "s/celldm[(]1[)]=[0-9]+/celldm(1)=${val}/"`

	elif [[ $param == 'k_point' ]]; then

	    
		 content=`echo "$content" |  sed -r "s/([0-9]+[[:blank:]][0-9]+[[:blank:]][0-9]+[[:blank:]])+/$val $val $val /"`
	    
	fi

	   
    done

    
    echo -e  "$content" > "tmp.sh"
   
    if [[ -d $struct ]]; then
	
	files=$(ls $struct)
	
    else

	files=$struct
	
    fi
    
   
    for struct in $files; do
	
	if [[ $2 == 'singleton' ]]; then
	    sbatch --time=$time --job-name=$Name --dependency=singleton "tmp.sh" $Name $struct
	elif [[ $option == 'epsil' ]]; then
	
	    sbatch --time=$time --job-name=$Name "tmp.sh" $Name $struct "kgrid.dat"
	
	else
	    sbatch --time=$time --job-name=$Name "tmp.sh" $Name $struct
	fi

    done
    
    rm  "tmp.sh"
}


Name=''
lat=0
struct=''
charge=0
time=0-00:05:00
ecut=0
k_point=0
options=0

declare -a modify
while getopts ':o:k:e:s:n:t:c:a:' c; do
    case $c in
	o) options=$OPTARG ;;
	k) k_point=$OPTARG; modify+=( 'k_point' );;
	e) ecut=$OPTARG;    modify+=( 'ecut' );;
	s) struct=$OPTARG;  modify+=( 'struct' )  ;;
	n) Name=$OPTARG; modify+=( 'Name' );           ;;
	t) time=$OPTARG;   modify+=( 'time' )        ;;
	c) charge=$OPTARG; modify+=( 'charge' )        ;;
	a) lat=$OPTARG; modify+=( 'lat' )         ;;
	-)
	    case ${OPTARG} in
		
		"help"|"Help") echo "not implemented" ;;
	    esac
			       
    esac

done


if [[ $Name == '' || $struct == '' ]]; then
    echo "Name and structure are madetory please provide using -n and -s"
    exit
fi

check_regex='[0-9][0-9]*\.*[0-9]*[-:,][0-9][0-9]*\.*[0-9]*[-:,][0-9][0-9]*\.*[0-9]*'
range_regex='s/\([0-9][0-9]*\.*[0-9]*\)[-:,]\([0-9][0-9]*\.*[0-9]*\)[-:,]\([0-9][0-9]*\.*[0-9]*\)/\1 \2 \3/'

name=$Name
if [[ $k_point =~ $check_regex ]]; then
    for val in $(seq `echo "$k_point" | sed "$range_regex"`); do

	k_point=$val
	Name="${name}_${k_point}"
	
	changeParams jobs/scf.sh
	sleep 2
	
    done
    
elif [[ $ecut =~ $check_regex ]]; then
    
    for val in $(seq `echo $ecut | sed "$range_regex"`); do
	ecut=$val
	Name="${name}_${ecut}"
	changeParams jobs/scf.sh
	sleep 2
	
    done
    
#requieres celldm(1) to be specified
elif [[ $lat =~ $check_regex ]]; then
    
    for val in $(seq `echo "$lat" | sed "$range_regex"`); do
	
	lat=$val
	Name="${name}_${lat}"
	
	changeParams jobs/scf.sh
	sleep 2
	
    done
    
elif [[ $options == 'bands' ]]; then

    
    changeParams jobs/scf.sh 'singleton'
    sleep 2

    changeParams jobs/band.sh 'singleton'
    sleep 2
    
    changeParams jobs/plot_bands.sh 'singleton'

elif [[ $options == 'relax' ]]; then
    
    changeParams jobs/relax.sh
    
elif [[ $options == 'irelax' ]]; then
    
    changeParams jobs/ion_relax.sh
    
fi
