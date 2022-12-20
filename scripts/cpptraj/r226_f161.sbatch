#!/bin/bash
 
#SBATCH --job-name=r226_f161
#SBATCH --output=r226_f161.out
#SBATCH --error=r226_f161.err
#SBATCH --partition=dinner
#SBATCH --account=pi-dinner
#SBATCH --qos=dinner
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6G

module load parallel
module load amber

ulimit -u 10000

srun="srun --exclusive -N1 -n1"

parallel="parallel --delay 0.2 -j $SLURM_NTASKS --joblog runtask.log --resume"

$parallel "$srun cpptraj -i rf.{1}.run" ::: {0..294}
