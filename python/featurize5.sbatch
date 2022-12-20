#!/bin/bash
#SBATCH --job-name=featurize
#SBATCH --output=featurize.out
#SBATCH --error=feat.err

#SBATCH --partition=rouxgpu
#SBATCH --account=pi-roux
#SBATCH --qos=roux

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=128G

#SBATCH --time=20:00:00

#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=$USER@uchicago.edu

source ~/.bashrc
zsh
conda activate py37

python featurize5.py
