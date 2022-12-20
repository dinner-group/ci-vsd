# run a VMD tcl script to save specificed frames from specified trajectories as pdb files

import numpy as np
import os
import sys
import subprocess

indices = np.loadtxt("../suppl-runs/seeds.txt", dtype=int)

for i, (traj, frame) in enumerate(indices):
    trajfile = f"../anton/xtc1000ns/civsd-{traj}.xtc"
    out = f"../suppl-runs/{i}.pdb"
    command = f"vmd -dispdev none -e save-frame.tcl -args {trajfile} {frame} {out}"
    # os.system(command)
    print(f"{i}: Extracting frame {frame} from trajectory {traj}")
    subprocess.run(command, shell=True, stdout=subprocess.DEVNULL)
