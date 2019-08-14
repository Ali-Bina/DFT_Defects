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

cat > ./inputs/nscf_${1}.in << EOF
&control
    calculation='nscf'
    restart_mode='from_scratch',
    prefix='GeTe'
    pseudo_dir='./PP'
    verbosity='high'
    tprnfor=.true.
    tstress=.true.
    outdir='./tmp'
/
&system
    ibrav=0,
    nat=2 , ntyp= 2, 
    celldm(1)=8.183,
    ecutwfc = 100, nbnd=50
    occupations='tetrahedra'
    nosym=.true.
    noncolin=.true.,
    lspinorb=.true.
/
&electrons
    diagonalization='david'
    conv_thr = 1.0e-10
    mixing_beta = 0.7
    electron_maxstep = 500
    mixing_mode = 'plain'
/
ATOMIC_SPECIES
Ge  72.630  Ge.rel-pbe-n-nc.UPF
Te  127.60  Te.rel-pbe-n-nc.UPF
CELL_PARAMETERS (alat=  8.18300000)
   0.489850122  -0.282815099   0.842668840
  -0.000000000   0.565630198   0.842668840
  -0.489850122  -0.282815099   0.842668840
ATOMIC_POSITIONS (crystal)
Ge       0.234772997   0.234772997   0.234772997
Te       0.765227003   0.765227003   0.765227003
K_POINTS (automatic)
30 30 30 0 0 0
EOF

srun pw.x < ./inputs/nscf_${1}.in > ./outputs/nscf_${1}.out