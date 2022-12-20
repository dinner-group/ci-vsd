#!/bin/bash

for i in $(seq -f "%03g" $1 $2); do
    cd $i/
    sbatch run.sbatch
    cd ../
done
