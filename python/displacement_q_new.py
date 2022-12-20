import numpy as np
import MDAnalysis as mda
import MDAnalysis.transformations
from joblib import Parallel, delayed
from multiprocessing import cpu_count

TOPFILE = "/project/dinner/scguo/ci-vsd/models/MD-clustering-center/civsd.psf"


def compute_displacement_charge(u, not_ions, ions):
    Lz = u.trajectory.dimensions[2]  # length of box
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


def traj_displacement_charge(topfile, trajfile, ion_sel="segid ION", n_frame=10000):
    u = mda.Universe(topfile, trajfile)

    # fix periodic boundary conditions
    all = u.atoms
    trans = mda.transformations.wrap(all, compound="residues")
    ions = u.select_atoms(ion_sel)
    # but not for ions
    ion_unwrap = mda.transformations.unwrap(ions)
    # center at (0, 0, 0) so that Rong's calculation works
    center = mda.transformations.center_in_box(all, point=(0, 0, 0), wrap=True)
    u.trajectory.add_transformations(trans, ion_unwrap, center)

    not_ions = u.select_atoms(f"not {ion_sel}")
    ions = u.select_atoms(ion_sel)
    q = np.zeros(n_frame)
    for i, _ in enumerate(u.trajectory[:n_frame]):
        q[i] = compute_displacement_charge(u, not_ions, ions)
    return q


def main():
    n_jobs = cpu_count()

    dcds = []
    for i in range(3, 119):
        if i == 82:
            continue
        dcds.append(f"/beagle3/dinner/scguo/anton2-backup/dcdfiles/civsd.{i}.dcd")
    print(len(dcds))

    ion_sel = "name SOD CLA"
    displacement_charges_2 = Parallel(n_jobs=n_jobs, verbose=20)(
        delayed(traj_displacement_charge)(TOPFILE, dcd, ion_sel=ion_sel) for dcd in dcds
    )
    print("Finished with Anton trajectories")
    np.save(
        "/project/dinner/scguo/ci-vsd/data/raw_feat/displacement_q_anton2.npy",
        displacement_charges_2,
    )


if __name__ == "__main__":
    main()
