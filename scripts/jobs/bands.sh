#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --time=0-02:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dftali0@gmail.com
#SBATCH --job-name=GeTe_relax
#SBATCH --mem=150G
#SBATCH --account=rrg-maassenj
#SBATCH --output=%x-%j.out

module purge
module load nixpkgs/16.09 gcc/7.3.0 openmpi/3.1.2 
module load quantumespresso/6.4

cat > ./inputs/bands_${1}.in << EOF
&BANDS
    prefix  = 'GeTe'
    outdir = './'
    filband = 'GeTe_bands.dat'
/
EOF

srun bands.x < ./inputs/bands_${1}.in >  ./outputs/bands_${1}.out
