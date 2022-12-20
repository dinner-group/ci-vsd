#!/bin/bash

for i in $(seq -f "%03g" $1 $2); do
    cd $i/
    cp -f ../run.sbatch .
    sbatch run.sbatch
    cd ../
done
