#!/bin/bash
#SBATCH --job-name=smd-prep
#SBATCH --out=smd.out
#SBATCH --error=smd.err

#SBATCH --partition=beagle3
#SBATCH --account=pi-dinner

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --gres=gpu:1
#SBATCH --mem=32G
#SBATCH --constraint=a100

#SBATCH --time=24:00:00

export PLUMED_KERNEL=/software/plumed-2.7-el8-x86_64/lib/libplumedKernel.so

# load python
source ~/.bashrc
source ~/.zshrc
conda activate py39

trajfile=$1
frame=$2
outfolder=$3

echo "Trajectory file: $trajfile"
echo "Frame: $frame"
echo "Output folder: $outfolder"

echo "Prepping input files..."
# make input files
python prep.py $trajfile $frame $outfolder
# copy over necesary files
cp run33/civsd.prmtop run33/smd.mdin $outfolder

cd $outfolder
echo "Trajectory file: $trajfile" >> out.log
echo "Frame: $frame" >> out.log
echo "Output folder: $outfolder" >> out.log

#  prepare PLUMED file
pdbfile=civsd.pdb
nowatfile=civsd-nowat.pdb
# target_pdb=/project/dinner/scguo/ci-vsd/biased/smd-prep/down-seed2-nowat-noh.pdb
target_pdb=/project/dinner/scguo/ci-vsd/biased/smd-prep/up-seed2-nowat-noh.pdb
# target_pdb=/project/dinner/scguo/ci-vsd/biased/smd-prep/downminus-seed-nowat-noh.pdb

echo "Target pdb: $target_pdb" >> out.log

# 1st/2nd argument are reference pdbs
# 3rd is mode, "compute" or "pull"
# 4th is PLUMED file
# 5th is output of PLUMED calculations
# 6th/7th are start/end CV files for SMD 
module unload openmpi amber
module load plumed 
python ../plumed.py $nowatfile $target_pdb compute plumed_0.dat cv-start.txt
plumed driver --mf_pdb $nowatfile --plumed plumed_0.dat
python ../plumed.py $nowatfile $target_pdb compute plumed_1.dat cv-end.txt
plumed driver --mf_pdb $target_pdb --plumed plumed_1.dat
python ../plumed.py $nowatfile $target_pdb pull plumed.dat cv-smd.txt cv-start.txt cv-end.txt

# run initial SMD
module unload plumed intelmpi
module load amber
echo "Running SMD simulation..."
echo "    Started at `date`"
cpptraj -p civsd.prmtop -y $pdbfile -x civsd.rst7
pmemd.cuda -O -i smd.mdin -p civsd.prmtop -c civsd.rst7 -o amber.log -r restart0.rst7 -x run0.nc 
echo "Finished SMD simulation."
echo "    Finished at `date`"

# convert restart file to pdb
echo "Preparing new files..."
cpptraj -p civsd.prmtop -x restart0.pdb -y restart0.rst7
# strip waters and zero occupancy for second run
rest_nowat=restart0-nowat.pdb
python ../strip_pdb.py restart0.pdb $rest_nowat

# repeat PLUMED file wrangling
module unload amber openmpi
module load plumed

echo "New target: $nowatfile"
echo "New starting file: $rest_nowat"
python ../plumed.py $rest_nowat $nowatfile compute plumed_2.dat cv-restart.txt
plumed driver --mf_pdb $rest_nowat --plumed plumed_2.dat
# new target is original PDB file
python ../plumed.py $rest_nowat $nowatfile compute plumed_3.dat cv-reend.txt
plumed driver --mf_pdb $nowatfile --plumed plumed_3.dat
python ../plumed.py $rest_nowat $nowatfile pull plumed.dat cv-resmd.txt cv-restart.txt cv-reend.txt

# run backwards pulling SMD
module unload plumed intelmpi
module load amber

echo "Running SMD simulation..."
echo "    Started at `date`"
pmemd.cuda -O -i smd.mdin -p civsd.prmtop -c restart0.rst7 -o amber.log -r restart1.rst7 -x run1.nc 
echo "Finished SMD simulation."
echo "    Finished at `date`"


