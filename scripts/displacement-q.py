import sys
import os
import glob

import numpy as np
import MDAnalysis as mda
from MDAnalysis import transformations

from joblib import Parallel, delayed
from multiprocess import cpu_count

n_jobs = cpu_count()


HOME_DIR = "/project/dinner/scguo/ci-vsd"
TOPFILE = f"{HOME_DIR}/models/MD-clustering-center/civsd.psf"


def compute_displacement_charge(u, not_ions, ions):
    Lz = u.dimensions[2]  # length of box
    q = 0
    # protein, lipid, water
    for atom in not_ions:
        q += atom.charge * (atom.position[2] + 0.5 * Lz) / Lz
        # ions, use unwrapped coordinates
    for atom in ions:
        z = atom.position[2]
        if z > 0:
            z -= Lz
        q += atom.charge * (z + 1.0 * Lz) / Lz
    return q


def compute_displacement_charge_single(u, not_ions, ions):
    Lz = u.dimensions[2]  # length of box
    q = 0
    # protein, lipid, water
    # for atom in u.select_atoms("not segid ION"):
    for atom in not_ions:
        q += atom.charge * (atom.position[2] + 0.5 * Lz) / Lz
    # ions, use unwrapped coordinates
    # for atom in u.select_atoms("segid ION"):
    for atom in ions:
        z = atom.position[2]
        if z > 0:
            z -= Lz
        q += atom.charge * (z + 1.0 * Lz) / Lz
    return q


def traj_displacement_charge(
    topfile, trajfile, ion_sel="segid ION", n_frame=10000
):
    u = mda.Universe(topfile, trajfile)

    # wrap and center
    prot = u.select_atoms("protein")
    trans = transformations.unwrap(prot)
    center = transformations.center_in_box(prot, point=(0, 0, 0), wrap=False)
    u.trajectory.add_transformations(trans, center)

    # compute dq
    not_ions = u.select_atoms(f"not {ion_sel}")
    ions = u.select_atoms(ion_sel)
    q = np.zeros(n_frame)
    for i, _ in enumerate(u.trajectory[:n_frame]):
        q[i] = compute_displacement_charge(u, not_ions, ions)
        # q[i] = compute_displacement_charge_nopbc(u, not_ions, ions)
    return q


def main():
    water_dcds = []
    results_first = np.load(f"{HOME_DIR}/data/raw_feat/displacement_q_anton2_start.npy")
    # results_first = np.zeros((115))
    for i in range(3, 119):
        if i == 82:
            continue
        water_dcds.append(f"/beagle3/dinner/scguo/anton2-backup/dcdfiles/civsd.{i}.dcd")

    result = Parallel(n_jobs=n_jobs, verbose=20)(
        delayed(traj_displacement_charge)(
            TOPFILE, dcd, ion_sel="name SOD CLA"
        ) for dcd in water_dcds
    )
    displacement_trajs = np.zeros((115, 10001))
    displacement_trajs[:, 1:] = result
    displacement_trajs[:, 0] = results_first
    np.save("../data/raw_feat/displacement_q_anton2.npy", displacement_trajs)


if __name__ == "__main__":
    main()
