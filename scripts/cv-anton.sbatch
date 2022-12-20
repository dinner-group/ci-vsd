#!/bin/bash
#SBATCH --job-name=vmd
#SBATCH --output=dist-spin.out
#SBATCH --error=vmd.err

#SBATCH --partition=roux 
#SBATCH --account=pi-roux
#SBATCH --qos=roux

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G

#SBATCH --time=08:00:00

#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=$USER@uchicago.edu


module load vmd
for i in {0..294};
do
    vmd -dispdev text -e dist-spin.tcl -args civsd-${i} 
done
