import sys
import os
import glob

import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis

from joblib import Parallel, delayed
from multiprocess import cpu_count

n_jobs = cpu_count()


HOME_DIR = "/project/dinner/scguo/ci-vsd"
CUTOFF = 3.5
ANGLE = 120


def hbond_anton(startfile, trajfile, r_i, cutoff=3.4, angle=90):
    u = mda.Universe(startfile)
    u = u.load_new(trajfile)
    r_sel = f"protein and resid {r_i} and not backbone and type N"
    phos_sel = f"name O11 O12 O13 O14 and around 3.0 (protein and resid {r_i})"
    h_sel = f"protein and resid {r_i} and name HE HH11 HH21 HH22 HH12"
    hbonds = HydrogenBondAnalysis(
        universe=u,
        donors_sel=r_sel,
        hydrogens_sel=h_sel,
        acceptors_sel=phos_sel,
        update_selections=True,
        d_a_cutoff=cutoff,
        d_h_a_angle_cutoff=angle,
    )
    return hbonds


def hbond_old(trajfile, r_i, cutoff=3.4, angle=90):
    u = mda.Universe(f"{HOME_DIR}/models/MD-clustering-center/civsd.psf", trajfile)
    r_sel = f"protein and resid {r_i} and not backbone and type NC2"
    phos_sel = f"name O11 O12 O13 O14 and around 3.0 (protein and resid {r_i})"
    h_sel = f"protein and resid {r_i} and type HC"
    hbonds = HydrogenBondAnalysis(
        universe=u,
        donors_sel=r_sel,
        hydrogens_sel=h_sel,
        acceptors_sel=phos_sel,
        update_selections=True,
        d_a_cutoff=cutoff,
        d_h_a_angle_cutoff=angle,
    )
    return hbonds


def run_analysis(analysis):
    analysis.run()
    return analysis.count_by_time()


def compute_hb_cutoff_angle(hb_fn, r_i, startfiles, trajfiles, cutoff, angle, n_jobs=40, verbose=10):
    analysis_ensemble = [hb_fn(s, t, r_i, cutoff=cutoff, angle=angle) for (s, t) in zip(startfiles, trajfiles)]
    results = Parallel(n_jobs=n_jobs, verbose=verbose)(
        delayed(run_analysis)(analysis) for analysis in analysis_ensemble
    )
    return results


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
    dcds = []
    startfiles = []
    for i in range(3, 119):
        if i == 82:
            continue
        dcds.append(f"{HOME_DIR}/anton2/nowat/civsd.nowat.{i}.dcd")
        dmsfiles = glob.glob(f"/beagle3/dinner/scguo/anton2-backup/workdir.{i}/*.dms")
        k = ""
        for f in dmsfiles:
            if "groups" not in f:
                k = f
        startfiles.append(k)

    # dcds_old = []
    # for i in range(179):
    #     dcds_old.append(f"/project/dinner/scguo/anton-old/civsd_{i}.dcd")

    # print(f"Computing hydrogen bonds with cutoff {CUTOFF} and angle {ANGLE}")
    results = []
    # results_old = []
    # dcds_second = files_second()
    # results_second = []
    for r_i in (217, 223, 226, 229, 232):
        results.append(
            compute_hb_cutoff_angle(
                hbond_anton, r_i, [None]*115, dcds, CUTOFF, ANGLE, n_jobs=n_jobs, verbose=10
            )
        )
        # single_old = compute_hb_cutoff_angle(
        #     hbond_old, r_i, dcds_old, CUTOFF, ANGLE, n_jobs=n_jobs, verbose=10
        # )
        # single_old = [t[:10000] for t in single_old]
        # single_second = compute_hb_cutoff_angle(
        #     hbond_old, r_i, dcds_second, CUTOFF, ANGLE, n_jobs=n_jobs, verbose=10
        # )
        # hb_new_short = []
        # for hb in single_second:
        #     length = len(hb)
        #     if length > 10000 and length < 100000:
        #         length = 10000
        #     hb_new_short.append(hb[:length])
        # results_second.append(hb_new_short)
    hb_phos = np.stack(results, axis=-1)
    print(hb_phos.shape)
    np.save(f"{HOME_DIR}/data/raw_feat/hbond_phos_anton2.npy", hb_phos)

    # hb_phos_old = np.stack(results_old, axis=-1)
    # np.save(f"{HOME_DIR}/data/raw_feat/hbond_phos_0-179.npy", hb_phos_old)

    # hb_phos_second = np.stack(results_second, axis=-1)
    # np.save(f"{HOME_DIR}/data/raw_feat/hbond_phos_179-300.npy", results_second)


if __name__ == "__main__":
    main()
