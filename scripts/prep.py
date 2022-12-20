import os
import subprocess
import numpy as np

data = np.loadtxt("seeds.txt", dtype=int)
for i, (traj, frame) in enumerate(data):
    newdir = str(i).zfill(3)
    if not os.path.exists(newdir):
        os.mkdir(newdir)
    subprocess.call(
        f"cp -f anton2-extract-pdbs/civsd-{traj}-{frame}.pdb {newdir}/civsd.pdb",
        shell=True,
    )
    subprocess.call(f"cp -f test/civsd.prmtop {newdir}", shell=True)
    subprocess.call(f"cp -f test/npt.mdin {newdir}", shell=True)
    subprocess.call(f"cp -f test/run.sbatch {newdir}", shell=True)
