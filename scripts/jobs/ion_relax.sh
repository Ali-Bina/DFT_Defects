#!/bin/bash

#SBATCH --nodes=2
#SBATCH --ntasks-per-node=40
#SBATCH --time=0-02:00:00
#SBATCH --mem=150G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dftali0@gmail.com
#SBATCH --job-name=GeTe_relax
#SBATCH --account=rrg-maassenj
#SBATCH --output=%x-%j.out

module purge
module load nixpkgs/16.09 gcc/7.3.0 openmpi/3.1.2 
module load quantumespresso/6.4


nat=$(( `cat $2 | wc -l` - 3 ))

cat > ./inputs/relax_${1}.in <<EOF

&control
    calculation='relax'
    restart_mode='from_scratch',
    prefix='$1'
    pseudo_dir='./PP'
    verbosity='high'
    tprnfor=.true.
    nstep=500
    tstress=.true.
    etot_conv_thr=1.0e-04
    forc_conv_thr=1.0e-03
    dt=20
    outdir='./tmp'	
/
&system
    ibrav=0,
    nat=$nat, ntyp= 2,
    ecutwfc=100,
    occupations='fixed',
    noncolin=.true.,
    lspinorb=.true.
    tot_charge=0	
    nbnd=372
/
&electrons
    diagonalization='david'
    conv_thr = 1.0e-6
    mixing_beta = 0.5
    electron_maxstep = 500
    mixing_mode = 'plain'
/
&IONS
  ion_dynamics = 'bfgs'
/
ATOMIC_SPECIES
Ge  72.630  Ge.rel-pbe-n-nc.UPF
Te  127.60  Te.rel-pbe-n-nc.UPF

CELL_PARAMETERS (bohr)
`head -n 3 $2`

ATOMIC_POSITIONS (crystal)
`tail -n $nat $2`

K_POINTS (automatic) 
4 4 5  0 0 0

EOF

srun pw.x < ./inputs/relax_${1}.in > ./outputs/relax_${1}.out
