#!/bin/bash
 
#SBATCH --job-name=hbonds
#SBATCH --output=hbonds.out
#SBATCH --error=hbonds.err
#SBATCH --partition=dinner
#SBATCH --account=pi-dinner
#SBATCH --qos=dinner
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=357G

cd /project/dinner/scguo
# load Anaconda here
source ~/.bashrc
zsh

env=py39
conda activate $env

# script="/project/dinner/scguo/ci-vsd/scripts/hbonds-lipid.py"
# script="/project/dinner/scguo/ci-vsd/scripts/hbonds-water.py"
script="/project/dinner/scguo/ci-vsd/scripts/hbonds-sb.py"

python $script
