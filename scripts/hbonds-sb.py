import sys
import os

import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis

from joblib import Parallel, delayed
from multiprocess import cpu_count

n_jobs = cpu_count()


HOME_DIR = "/project/dinner/scguo/ci-vsd"
TOPFILE = f"{HOME_DIR}/civsd-pro.psf"
CUTOFF = 3.5
ANGLE = 120


def hbond_traj(trajfile, cutoff=3.4, angle=90):
    u = mda.Universe(TOPFILE, trajfile)
    r_sel = f"protein and resid 217 223 226 229 232 and name NE NH1 NH2"
    d_sel = f"protein and resid 129 183 186 and name OD1 OD2 OE1 OE2"
    hbonds = HydrogenBondAnalysis(
        universe=u,
        donors_sel=r_sel,
        hydrogens_sel=None,
        acceptors_sel=d_sel,
        update_selections=False,
        d_a_cutoff=cutoff,
        d_h_a_angle_cutoff=angle,
    )
    return hbonds


def run_analysis(analysis):
    analysis.run()
    return analysis


def list_files():
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
            files.append(f"/project/dinner/scguo/anton-old/xtc1000ns/civsd-{i}.xtc")
    for i in range(3, 119):
        if i == 82:
            continue
        files.append(f"{HOME_DIR}/anton2/prot/civsd.prot.{i}.xtc")
    return files


def main():
    dcds = list_files()
    all_analyses = [hbond_traj(f, cutoff=CUTOFF, angle=ANGLE) for f in dcds]
    results_all = Parallel(n_jobs=n_jobs, verbose=20)(
        delayed(run_analysis)(analysis) for analysis in all_analyses)

    # get index of ids
    u = mda.Universe(TOPFILE, dcds[0])
    arrs = []
    # array goes R217-D129, R217-D183, R217-D186, R223-D129, ..., R232-D186
    index_dic = {}
    for i, r in enumerate((217, 223, 226, 229, 232)):
        for j, d in enumerate((129, 183, 186)):
            index_dic[(r, d)] = i * 3 + j
    for i in results_all:
        result = i.results.hbonds
        t_len = len(i.frames)
        arr = np.zeros((t_len, 18), dtype=int)
        for (frame, d_id, _, a_id, _, _) in result:
            d_resid = u.atoms[int(d_id)].resid
            a_resid = u.atoms[int(a_id)].resid
            y = index_dic[(d_resid, a_resid)]
            arr[int(frame), y] += 1
        arrs.append(arr)
    np.save(f"{HOME_DIR}/data/raw_feat/hbonds_sb.npy", arrs[:237])
    np.save(f"{HOME_DIR}/data/raw_feat/hbonds_sb_anton2.npy", arrs[237:])


if __name__ == "__main__":
    main()
