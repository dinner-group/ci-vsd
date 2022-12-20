import mdtraj as md
import glob
import os
import sys

start = int(sys.argv[1])
end = int(sys.argv[2])
for i in range(start, end + 1):
    print(f"Step {i}")
    os.chdir(f"./{str(i).zfill(3)}")
    trjin = sys.argv[3]
    trjout = sys.argv[4]
    trj = md.load_netcdf(trjin, top="civsd.prmtop")
    trj[-3000:].atom_slice(trj.topology.select("protein")).save_xtc(trjout)
    os.chdir("..")

# for i in range(150, 200):
#    os.chdir('/project/roux/ymeng1/kinases/eli_lilly/abl_training/md/cluster%s' % str(i))
#    for j in range(13, 14):
#        ifile = 'abl_' + str(i) + '_' + str(j) + '.trj'
#        ofile = 'abl_' + str(i) + '_' + str(j) + '.xtc'
#        trj = md.load_netcdf(ifile, top='abl_%s.prmtop' % str(i))
#        trj.atom_slice(trj.topology.select('protein')).save_xtc(ofile)
