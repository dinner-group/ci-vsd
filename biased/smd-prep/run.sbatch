#!/bin/bash
#SBATCH --job-name=smd-amber
#SBATCH --out=smd.out
#SBATCH --error=smd.err

#SBATCH --partition=gm4
#SBATCH --account=pi-dinner
#SBATCH --qos=gm4

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1

#SBATCH --time=7:00:00

#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=$USER@uchicago.edu

module purge
module load amber
module unload flex
module load plumed/2.3

export PLUMED_KERNEL=/software/plumed-2.3-el7-x86_64/lib/libplumedKernel.so

cpptraj -p civsd.prmtop -y civsd.pdb -x civsd.rst7
pmemd.cuda -O -i smd.mdin -p civsd.prmtop -c civsd.rst7 -o amber.log -r civsd.rst7 -x run0.nc 
