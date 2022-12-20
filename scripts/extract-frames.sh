#!/bin/bash

# load VMD
module load vmd

i=0
while IFS= read -r line; do
    IFS=' ' read -r -a array <<< "$line"
    traj=${array[0]}
    frame=${array[1]}
    outfile=../suppl-runs/${i}.pdb
    trajfile=../anton/xtc1000ns/civsd-${traj}.xtc
    
    vmd -dispdev none -e save-frame.tcl -args $trajfile $frame $outfile

    i=$((i+1))
done < ../suppl-runs/seeds.txt
