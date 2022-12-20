#!/bin/bash

# must have python loaded first
# $1 is first trajectory number
# $2 is 2nd trajectory number
# $3 is input trajectory file name
# $4 is water-free trajectory output file name
echo "Extracting protein..."
python extract_protein.py $1 $2 $3 $4

# load vmd
module load vmd/1.9.3

# $5 is output CV file
echo "Calculating CVs..."
for i in $(seq -f "%03g" $1 $2); do
    vmd -dispdev none -e dist-spin.tcl -args ${i}/$4 ${i}/$5
done
