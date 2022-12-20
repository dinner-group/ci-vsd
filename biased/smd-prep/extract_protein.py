import mdtraj as md
import glob
import os
import sys

infile = sys.argv[1]
outfile = sys.argv[2]
print(f"Input trajectory file {infile}")
print(f"Output trajectory file {outfile}")

trj = md.load_netcdf(infile, top="civsd.prmtop")
trj.atom_slice(trj.topology.select("protein")).save_xtc(outfile, force_overwrite=True)
