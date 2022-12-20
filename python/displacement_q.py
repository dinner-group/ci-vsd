import numpy as np
import MDAnalysis as mda
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


def compute_displacement_charge_nopbc(u, not_ions, ions):
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
        q += atom.charge * (z + 1.0 * Lz) / Lz
    return q


def traj_displacement_charge(topfile, trajfile, ion_sel="segid ION", n_frame=10000):
    u = mda.Universe(topfile, trajfile)
    not_ions = u.select_atoms(f"not {ion_sel}")
    ions = u.select_atoms(ion_sel)
    q = np.zeros(n_frame)
    for i, _ in enumerate(u.trajectory[:n_frame]):
        q[i] = compute_displacement_charge(u, not_ions, ions)
        # q[i] = compute_displacement_charge_nopbc(u, not_ions, ions)
    return q


def files_second():
    remove = {
        1282,
        1283,
        1284,
        1285,
        1286,
        1288,
        1289,
        1290,
        1187,
        1188,
        1189,
        1190,
        1191,
        1197,
        1198,
        1199,
        1203,
        1205,
        1206,
        1207,
        1211,
        1212,
        1213,
        1214,
        1215,
        1225,
        1226,
        1227,
        1228,
        1231,
        1232,
        1233,
        1236,
        1237,
        1238,
        1242,
        1245,
        1246,
        1252,
        1253,
        1260,
        1261,
        1262,
        1263,
        1266,
        1267,
        1268,
        1269,
        1270,
        1271,
        1272,
        1273,
        1274,
        1275,
        1276,
        1277,
        1278,
        1279,
    }
    files = []
    for i in range(179, 295):
        if i == 180:
            continue
        if (i + 1000) not in remove:
            files.append(f"/project/dinner/scguo/anton-old/civsd_{i}.dcd")
    return files


def main():
    n_jobs = cpu_count()

    # dcds_old = []
    # for i in range(179):
    #     dcds_old.append(f"/project/dinner/scguo/anton-old/civsd_{i}.dcd")
    # print(len(dcds_old))

    # displacement_charges = Parallel(n_jobs=n_jobs, verbose=20)(
    #     delayed(traj_displacement_charge)(TOPFILE, dcd) for dcd in dcds_old)
    # print(f"Finished with old Anton trajectories 0-178")
    # np.save("/project/dinner/scguo/ci-vsd/data/raw_feat/displacement_q_0-178.npy", displacement_charges)

    # dcds = []
    # for i in range(3, 119):
    #     if i == 82:
    #         continue
    #     dcds.append(f"/beagle3/dinner/scguo/anton2-backup/dcdfiles/civsd.{i}.dcd")
    # print(len(dcds))

    # ion_sel = "name SOD CLA"
    # displacement_charges_2 = Parallel(n_jobs=n_jobs, verbose=20)(
    #     delayed(traj_displacement_charge)(TOPFILE, dcd, ion_sel=ion_sel) for dcd in dcds
    # )
    # print("Finished with Anton trajectories")
    # np.save(
    #     "/project/dinner/scguo/ci-vsd/data/raw_feat/displacement_q_anton2.npy",
    #     displacement_charges_2,
    # )

    dcds_second = files_second()
    displacement_charges = Parallel(n_jobs=8, verbose=20)(
        delayed(traj_displacement_charge)(TOPFILE, dcd, n_frame=100000) for dcd in dcds_second
    )
    disp_short = []
    for dq in displacement_charges:
        length = len(dq)
        if length > 10000 and length < 100000:
            length = 10000
        disp_short.append(dq[:length])
    np.save(
        "/project/dinner/scguo/ci-vsd/data/raw_feat/displacement_q_179-end.npy",
        disp_short,
    )


if __name__ == "__main__":
    main()
