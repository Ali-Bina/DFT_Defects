#!/bin/bash

##########################################################################################
#Script for submiting jobs to the cluster. Modifies parameters of a pre-exisiting input file (e.g. ecut, k-point, structure,...)
#Depending on the provided option (-o) perfroms a single scf calculation, band structure, relaxation, DOS, dielectric tensor or convergence w.r.t. Kpoint, lattice parameter or cut off energy
#Generates input files stored in the inputs director
#QE outputs are stored in the output directory
#PP directory stroes psuedos
#structures stores the crystal structure:
#        line 1  : lattice parameter
#        line 2-4: primitive lattice vectors
#        line 5- : atomic basis
# To perform a convergence wrt k, ecut or lattice param provide input in the format:
# start:step:stop
#This will submit a series of scf calculations by varying the value  for k/ecut/latParam in the specified range
##########################################################################################


#generate a uniform grid of kpoints in coordinates of the reciprocal lattice and modify QE input
kpoints(){
    inputFile=$4
    submit=''
    stop=`echo "$inputFile" | grep -n 'K_POINTS' | cut -d':' -f 1`
    while IFS= read -r line; do
	let i++
	submit=`echo -e  "$submit\n$line"`
	
	if (( $i >= $(( $stop - 1 ))  )); then
       
	    break
	
	fi 
	
	
    
    done <<< "$inputFile"
    nPoints=$(( $1 * $2 * $3  ))
    submit=`echo -e "$submit\nK_POINTS crystal\n$nPoints"`
    
    #equally weighted k grid in crstal coordinates
    for (( i=0; i < $1; i++)){
	for (( j=0; j < $2; j++)){
		for (( k=0; k < $3; k++)){
		        
			Point=` awk -v i=$i -v j=$j -v k=$k -v l=$1 -v m=$2 -v n=$3 'BEGIN {printf "%0.15f\t%0.15f\t%0.15f\t%0.15f\n", i/l, j/m, k/n, 1 / (n * m *l)}'`
			submit=`echo -e "$submit\n$Point"`
			
		    }
	    }
	}

     echo -e "$submit"
    
}


#Replaces old structure with new one
modifyStruct(){

    inputFile=$1 #contents of input file
    newStruct=$2 #File containing the new structure
    submit=''
    #copy over all Parameters from input up to ATOMIC_POSITIONS 
    i=0
    stop=`echo "$inputFile" | grep -n 'CELL_PARAMETERS' | cut -d':' -f 1`
    
    while IFS= read -r line; do
	let i++
	submit=`echo -e  "$submit\n$line"`
	
	if (( $i >= $stop )); then
       
	    break
	
	fi 
	
	
    
    done <<< "$inputFile"
    
    fileLen=`cat $newStruct | wc -l`
    nAtoms=$(( $fileLen  - 4  ))
    newCellParam=`head -n 4 $newStruct | tail -n 3`
    newCoords=`tail -n $nAtoms $newStruct`
    
    submit=`echo -e "$submit\n$newCellParam"`
    
    atCard=`echo -e  "$inputFile" | grep 'ATOMIC_POSITIONS'`
    submit=`echo -e "$submit\n$atCard"`
    submit=`echo -e "$submit\n$newCoords"`
    
    stop=`echo -e "$inputFile" | grep -n '^[[:blank:]]*[a-zA-Z][a-zA-Z]*[[:blank:]]*[0-9][0-9]*\..*[0-9][[:blank:]]*$' | cut -d':' -f 1 | tail -n 1`
    fileLen=`echo -e "$inputFile" | wc -l`
    remainder=`echo -e "$inputFile" | tail -n $(( $fileLen - $stop ))`
    submit=`echo -e "$submit\n$remainder"`
    submit=`echo -e "$submit" | sed -r "s/nat[[:blank:]]*=[[:blank:]]*[0-9][0-9]*/nat=$nAtoms/"`
    echo -e "$submit"
   
}

#modifies Parameters of interest from the input
changeParams(){
    content=$1
    for param in ${modify[@]}; do

	val=`eval echo \\${$param}`
	
	if [[ $param == 'charge'  ]]; then
	    
	    content=`echo "$content" | sed -r  "s/tot_charge[[:blank:]]*=[[:blank:]]*[-]*[0-9]+\.*[0-9]*/tot_charge=${val}/"`
	    
	elif [[ $param == 'ecut' ]];then

	    content=`echo "$content" | sed -r  "s/ecutwfc[[:blank:]]*=[[:blank:]]*[0-9]+\.*[0-9]*/ecutwfc=${val}/"`

	elif [[ $param == 'lat' ]]; then

	    content=`echo "$content" | sed -r  "s/celldm[(]1[)][[:blank:]]*=[[:blank:]]*[0-9]+\.*[0-9]*/celldm(1)=${val}/"`

	elif [[ $param == 'k_point' ]]; then

	    
	    content=`echo "$content" |  sed -r "s/([0-9]+[[:blank:]][0-9]+[[:blank:]][0-9]+[[:blank:]])+/$val $val $val /"`

	elif [[ $param == 'smearParam' ]]; then

	    content=`echo "$content | sed -r  "s/degauss[[:blank:]]*=[[:blank:]]*[0-9]+\.*[0-9]*/degauss=${val}/"`
	    
	fi

	   
    done

    echo -e "$content"
    
   
}

#concatenates scheduler directives with specified QE program to create a job script
makeJobScript(){
    #$1 = directive file, $2 = QE program $3 = input file name, $4 = output file name

    jobscript='#!/bin/bash'
    input=`cat $1`
    run="srun $2 < inputs/$3  > outputs/$4"
    jobscript=`echo -e "$jobscript\n$input\n$run"`    

    echo -e "$jobscript"
}


#modifies an input file, generates a temporary job script and submits it to the scheduler, writes
#new input file to inputs directory
submit(){
    #$1: input file to be modified
    #$2: QE program name: pw.x, ph.x,....
    #$3: if set to singleton, runs sbatch with --denpendcy=singleton (run jobs with the same name one after the other)
    
    input=`cat $1` #QE input file to be modified
    calc=$2 #calculation type

    #if a directory containing several structure files is specified
    if [[ -d $struct ]]; then

    	files=$(ls ${struct}/*)
	
    else

    	files=$struct
	
    fi
    
   #if a directory containing structure files is specified, the specified job is executed for each structure in that directory
    for struct in $files; do

	#modify structure if a structure file has been provided and the calculation uses pw.x
	if [[ $struct != "DNE" ]]; then
	    [ $calc == pw.x  ] &&  input=`modifyStruct "$input" "$struct"`
	fi
	
	input=`changeParams "$input"`

	
	if (( l != 0 && m != 0 && n != 0 )); then
	    
	    [ $calc == pw.x  ] && input=`kpoints $l $m $n "$input"`

	    #can't use symmetry with epsilon.x
	    check=`echo -e "$input" | grep 'nosym'`
	    if [ ! $check ]; then
		
		input=`echo  -e "$input" | sed -e 's/\&SYSTEM/\&SYSTEM\nnosym=.true./I'`

	    fi
	    
	fi
	
	inputName="${name}.in"
	echo -e "$input" > "inputs/$inputName" #write new input to inputs directory
	
	outputName="${name}.out"

	#Write temporary job script to current directory
	jobScript=`makeJobScript inputs/sch_dir $calc $inputName $outputName`
	echo -e "$jobScript" > "tmp.sh"

	##########testing code#########
	cat tmp.sh
	echo -e "$input"
	
	########testing code########
	
	#submit temporary jobscript to scheduler
	if [[ $3 == 'singleton' ]]; then
	 
    	    sbatch --time=$time --job-name=$Name --dependency=singleton "tmp.sh"
	    
    	else
	 
    	    sbatch --time=$time --job-name=$Name "tmp.sh"
    	fi

    done
    
    rm  "tmp.sh" #delete jobscript
}

k_point=0 #K grid density
n=0 #grid points along x
m=0 #grid points along y
l=0 #grid points along z

Name='' #QE prefix
name='' #File name for I/O
lat=0 #lattice parameter
struct='DNE' #structure file
charge=0 #charge state
time=0-00:05:00 
ecut=0 #kinetic energy cutoff
options=0 #type of calculation
smearing=0

declare -a modify #stores the parameters to be modified
#take input from the command line
while getopts ':o:k:e:s:n:t:c:a:-:m:' c; do
    case $c in
	o) options=$OPTARG ;;
	k) k_point=$OPTARG; modify+=( 'k_point' );;
	e) ecut=$OPTARG;    modify+=( 'ecut' );;
	s) struct=$OPTARG;  modify+=( 'struct' )  ;;
	n) Name=$OPTARG; modify+=( 'Name' );           ;;
	t) time=$OPTARG;   modify+=( 'time' )        ;;
	c) charge=$OPTARG; modify+=( 'charge' )        ;;
	a) lat=$OPTARG; modify+=( 'lat' )         ;;
	m) k_grid=$OPTARG;;
	-)
	    case ${OPTARG} in
		smearing*) smearing=`echo $OPTARG | cut -d' ' -f 3`; modify+=( 'smearing' ); echo "$OPTARG";; 
		"help") 
		
		cat << EOF
parameters: -o: Options for calculations
	    -k: kgrid density
    	    -e: Energy cutoff for basis
	    -a: lattice parameter
	    -s: structure file
	    -n: name of file
	    -t: execution time
      --manual: Specify a nxmxl kgrid manually (Rather than using QE's Automatic))
	--skip: skips the scf step for the calculations that require it
EOF
		exit;;
	    esac
	    
    esac
    
done


if [[ $Name == '' ]]; then
    echo -e "Name is madetory please provide using -n\nUse submit.sh --help for more information"
    exit
fi

#create required directories
if [[ ! -d outputs ]];then
    echo "Created a new input directory"
    mkdir outputs
fi

if [[ ! -e inputs ]];then
    echo "Created a new output directory"
    mkdir inputs
fi


check_regex='[0-9][0-9]*\.*[0-9]*[-:,][0-9][0-9]*\.*[0-9]*[-:,][0-9][0-9]*\.*[0-9]*'
range_regex='s/\([0-9][0-9]*\.*[0-9]*\)[-:,]\([0-9][0-9]*\.*[0-9]*\)[-:,]\([0-9][0-9]*\.*[0-9]*\)/\1 \2 \3/'

#extracts the number grid points along each direction, specified for a manual k-grid
declare -a gridPoints
if [[ $k_grid =~ $check_regex ]]; then
    
    for val in $(echo "$k_grid" | sed "$range_regex"); do
	gridPoints+=($val)
    done

    n=${gridPoints[0]}
    l=${gridPoints[1]}
    m=${gridPoints[2]}

fi



#Kpoint convergence
if [[ $k_point =~ $check_regex ]]; then
    for val in $(seq `echo "$k_point" | sed "$range_regex"`); do

	k_point=$val
	name="scf_${Name}_${k_point}"
	
	submit inputs/scf.in pw.x 
	sleep 2
	
    done

#ecut converge
elif [[ $ecut =~ $check_regex ]]; then
    
    for val in $(seq `echo $ecut | sed "$range_regex"`); do
	ecut=$val
	name="scf_${Name}_${ecut}"
	
	submit inputs/scf.in pw.x
	sleep 2
	
    done
#Smearing convergence
elif [[ $smearing =~ $check_regex ]]; then

    check=`grep 'degauss'`
    if [! $check ]; then 
	echo "Error degauss not specified in input"
	exit
    fi
    
    
      for val in $(seq `echo "$smearing" | sed "$range_regex"`); do

	smearing=$val
	name="scf_${Name}_${smearing}"
	
	submit inputs/scf.in pw.x 
	sleep 2
	
    done
    
#scf vs lattice Parameter requieres celldm(1) to be specified
elif [[ $lat =~ $check_regex ]]; then
    check=`grep 'celldm(1)'`
    if [! $check ]; then 
	echo "Error celldm(1) not specified in input"
	exit
    fi
    
    for val in $(seq `echo "$lat" | sed "$range_regex"`); do
	
	lat=$val
	name="scf_${Name}_${lat}"
	
	submit inputs/scf.in pw.x
	sleep 2
	
    done

#submits the necessary jobs required for a band structure calculation
elif [[ $options == 'bands' ]]; then

    name="scf_${Name}"
    submit inputs/scf.in pw.x 'singleton'
    sleep 2

    name="band_${Name}"
    submit inputs/band.in pw.x  'singleton'
    sleep 2

    name="bands_${Name}"
    submit inputs/bands.in bands.x 'singleton'

#submit a variable cell  relaxation
elif [[ $options == 'relax' ]]; then

    name="relax_${Name}"
    submit inputs/relax.in pw.x

#submits ionic relaxation
elif [[ $options == 'irelax' ]]; then

    name="relax_${Name}"
    submit jobs/ion_relax.in pw.x

#dielectric tensor
elif [[ $options == 'epsil' ]]; then

    #if a uniform equally weighted k grid is not provided ph.x is used to calculate
    if (( $n == 0 && $m == 0 && $l == 0 )); then

	echo "error: must specify a uniform k-grid using the -m option"

    fi

    name="scf_eps_${Name}"
    submit inputs/scf.in pw.x 'singleton'
    name="epsil_${Name}"
    submit inputs/epsil.in epsilon.x 'singleton'
    
fi
