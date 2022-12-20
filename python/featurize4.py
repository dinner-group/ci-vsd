# Featurize all trajectories suing feature 4 (indicators for
# contacts between S4 and S123)

import numpy as np
import pyemma
import glob


def main():
    # collect trajectories
    xtcs = []
    for i in range(0, 1000):
        xtcs.append(f"../amber/xtc300ns/civsd_{i}_300ns.xtc")
    for i in range(0, 295):
        xtcs.append(f"../anton/xtc1000ns/civsd-{i}.xtc")
    print(f"Number of trajectories to be loaded: {len(xtcs)}")

    # define features in pyemma
    feat4 = pyemma.coordinates.featurizer("../civsd-pro.pdb")
    pair_indices = np.loadtxt(
        "../amber-gpu/pyemma-civsd-feat5-msmb-combined/AtomIndices.txt", dtype=np.int32
    )
    feat4.add_residue_mindist(
        residue_pairs=pair_indices,
        scheme="closest-heavy",
        threshold=0.45,
        periodic=False,
    )
    print(f"Number of features: {len(feat4.describe())}")

    # load coordinates
    feat4_raw = pyemma.coordinates.load(xtcs, features=feat4, chunksize=64)

    # save coordinates
    np.savez_compressed("../data/raw_feat/feat4_raw.npz", feat4_raw)


if __name__ == "__main__":
    main()
