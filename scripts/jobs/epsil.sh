#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --time=00:15:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dftali0@gmail.com
#SBATCH --job-name=Si_nscf
#SBATCH --mem=150G
#SBATCH --account=rrg-maassenj
#SBATCH --output=%x-%j.out

module purge
module load nixpkgs/16.09 gcc/7.3.0 openmpi/3.1.2 
module load quantumespresso/6.4

cat > ./inputs/epsil.in <<EOF

&inputpp
outdir='./',
prefix='$1',
calculation='eps'
/

&energy_grid
smeartype='gauss',
intersmear=0.136,
intrasmear=0.0,
wmax=30.0,
wmin=0.0,
nw=600,
shift=0
/
EOF

srun epsilon.x < ./inputs/epsil.in > ./outputs/epsil.out
