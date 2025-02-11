#!/bin/bash

#SBATCH --nodes=2
#SBATCH --ntasks-per-node=40
#SBATCH --time=1-12:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dftali0@gmail.com
#SBATCH --job-name=scf3GeTe
#SBATCH --mem=150G
#SBATCH --account=rrg-maassenj
#SBATCH --output=%x-%j.out

module purge
module load nixpkgs/16.09 gcc/7.3.0 openmpi/3.1.2 
module load quantumespresso/6.4


nat=$(( `cat $2 | wc -l` - 4 ))
cat > ./inputs/scf_${1}.in << EOF

&control
    calculation='scf'
    restart_mode='from_scratch',
    prefix='$1'
    pseudo_dir='./PP'
    verbosity='high'
    tprnfor=.true.
    tstress=.true.
/
&system
    ibrav=0
    nat=$nat, ntyp= 2,
    ecutwfc=100,nbnd=370
    occupations = 'fixed' ,
    noncolin=.true.,
    lspinorb=.true.
/
&electrons
    diagonalization='cg'
    conv_thr = 1.0e-6
    mixing_beta = 0.5
    electron_maxstep = 500
    mixing_mode = 'plain'
/
ATOMIC_SPECIES
Ge  72.630  Ge.rel-pbe-n-nc.UPF
Te  127.60  Te.rel-pbe-n-nc.UPF
CELL_PARAMETERS (bohr)
`head -n 4 $2 | tail -n 3`

ATOMIC_POSITIONS (crystal)
`tail -n $nat $2`
K_POINTS {automatic} 
4 4 5  0 0 0
EOF
srun pw.x < ./inputs/scf_${1}.in > ./outputs/scf_${1}.out
