#!/bin/bash

mkdir $1
cp ./civsd.prmtop $1/
cp ./run.sbatch $1/
cp ./template_proj.dat $1/plumed.dat
cp ./smd.mdin $1/
cp ./down-seed2.pdb $1/civsd.pdb

cd $1/
sed -i -e "s/XX/$2/g" plumed.dat
sed -i -e "s/YY/$3/g" plumed.dat
sed -i -e "s/ZZ/$4/g" plumed.dat
sed -i -e "s/XXX/plumed.dat/g" smd.mdin

sbatch run.sbatch

cd ../
