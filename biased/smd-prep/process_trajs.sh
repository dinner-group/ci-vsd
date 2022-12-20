#!/bin/bash

# must have python loaded first
# $1 is input trajectory file name
# $2 is water-free trajectory output file name
echo "Extracting protein..."
python extract_protein.py $1 $2

# load vmd
module load vmd/1.9.3

# $3 is output CV file
echo "Calculating CVs..."
vmd -dispdev none -e dist-spin.tcl -args $2 $3
