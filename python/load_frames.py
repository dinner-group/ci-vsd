import os
import sys
import glob

import numpy as np
import mdtraj as md


def main(file_names, indices, top_file, outfile):
    # assemble file list
    xyz, angles, lengths = [], [], []
    for (f, idx) in zip(file_names, indices):
        frame = md.load_frame(f, idx, top=top_file)
        # if mdtraj encounters an impossible frame index it just
        # returns a 0-length trajectory
        if len(frame) != 0:
            xyz.append(frame.xyz)
            angles.append(frame.unitcell_angles)
            lengths.append(frame.unitcell_lengths)

    top = md.load_topology(top_file)
    traj = md.Trajectory(
        np.concatenate(xyz),
        top,
        unitcell_angles=np.concatenate(angles),
        unitcell_lengths=np.concatenate(lengths),
    )

    traj.save(outfile)


if __name__ == "__main__":
    index_file = sys.argv[1]
    top_file = sys.argv[2]
    outfile = sys.argv[3]

    if not os.path.exists(index_file):
        raise IOError(f"{index_file} does not exist")
    if not os.path.exists(top_file):
        raise IOError(f"{top_file} does not exist")

    print(f"Index file: {index_file}")
    print(f"Topology file: {top_file}")
    print(f"Output file: {outfile}")

    file_names = []
    indices = []
    with open(index_file, mode="r") as f:
        for line in f.readlines():
            name, ix = line.strip("\n").split()
            file_names.append(name)
            indices.append(int(ix))

    main(file_names, indices, top_file, outfile)
