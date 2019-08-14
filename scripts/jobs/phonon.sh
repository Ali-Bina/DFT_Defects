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


cat > ./inputs/phonon.in <<EOF

&inputph
outdir='./'
prefix='$1'
tr2_ph=1d-14
epsil=.true.
ldisp=.false.
amass(1)=72.630
amass(2)=127.60
fildyn='GeTe.dyn'
/

0 0 0

EOF
 


srun ph.x < ./inputs/phonon.in > ./outputs/phonon.out
