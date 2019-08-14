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

nat=$(( `cat $2 | wc -l` - 4 ))

cat > ./inputs/band_${1}.in << EOF

&control
    calculation     = 'bands'
    prefix          = 'GeTe'
    restart_mode    = 'from_scratch'
    wf_collect      = .false.
    pseudo_dir      = './PP'
    outdir          = './'
    verbosity       = 'high'
    tprnfor         = .true.
    tstress         = .true.
 /
 &system
    ibrav=0,
    nat=$nat , ntyp= 2, 
    celldm(1)=`head -n 1 $2`
    ecutwfc = 100, nbnd=50
    occupations='tetrahedra'
    noncolin=.true.,
    lspinorb=.true.
 /
 &electrons
    diagonalization = 'david'
    mixing_beta     = 0.3
    conv_thr        = 1.0d-9
 /
ATOMIC_SPECIES   
Ge  72.630  Ge.rel-pbe-n-nc.UPF
Te  127.60  Te.rel-pbe-n-nc.UPF
CELL_PARAMETERS (alat=`head -n 1 $2`)
 `head -n 4 $2 | tail -n 3`

ATOMIC_POSITIONS (crystal)
`tail -n $nat $2`

K_POINTS crystal_b
10
0.242 0.758 0.5  200 # W
0.5  0.5    0.5  200 # T
0.0 0.0 0.0 200 ! G
0.0 0.629 0.371 200 ! K
0.0  0.5  0.5  200 # X
0.0 0.0 0.0 200 ! G
0.0 0.5 0.0 200 ! L
0.371 0.758 0.371 200 ! U
0.2118 0.3714 0.1918 200 ! PI
0.0 0.0 0.0 1 ! G

EOF

srun pw.x < ./inputs/band_${1}.in > ./outputs/band_${1}.out
