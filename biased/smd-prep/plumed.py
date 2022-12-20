import sys
import os


def main(ref, mode, plumed, output, cv_files=None):
    """
    ref : list of 2 str
        Reference PDB structures
    mode : str "compute" or "pull"
        Whether to write CV calculation file or SMD pulling file
    plumed : str
    output : str
        Name of PLUMED output file
    cv_file : optional, list of 2 str
        Name of CV output files used to get target CV values for pulling.
        The first is for the starting structure, and the second for the
        target final structure.
    """
    ref1, ref2 = ref
    if mode == "compute":
        CV_STRING = f"""
WHOLEMOLECULES ENTITY0=1-2191

# compute RMSD to up+ state and down state for pushing/pulling
rmsd1: RMSD REFERENCE={ref1} TYPE=OPTIMAL
rmsd2: RMSD REFERENCE={ref2} TYPE=OPTIMAL

# salt bridges from ABMD analysis
sb0: DISTANCE ATOMS=1782,407  NOPBC # R217 - D129
sb1: DISTANCE ATOMS=1886,1300 NOPBC # R223 - D186
sb2: DISTANCE ATOMS=1948,531  NOPBC # R226 - D136
sb3: DISTANCE ATOMS=2004,407  NOPBC # R229 - D129
sb4: DISTANCE ATOMS=2057,1300 NOPBC # R232 - D186

PRINT ARG=* FILE={output} STRIDE=2500"""
        with open(plumed, mode="w+") as f:
            f.write(CV_STRING)
    elif mode == "pull":
        with open(cv_files[0], mode="r+") as f:
            # header line
            _ = next(f)
            lines = next(f)
        rmsd1s, rmsd2s, sb0s, sb1s, sb2s, sb3s, sb4s = lines.split()[1:8]
        with open(cv_files[1], mode="r+") as f:
            _ = next(f)
            linee = next(f)
        rmsd1e, rmsd2e, sb0e, sb1e, sb2e, sb3e, sb4e = linee.split()[1:8]

        SMD_STRING = f"""
WHOLEMOLECULES ENTITY0=1-2191

# compute RMSD to up+ state and down state for pushing/pulling
rmsd1: RMSD REFERENCE={ref1} TYPE=OPTIMAL
rmsd2: RMSD REFERENCE={ref2} TYPE=OPTIMAL

# salt bridges from ABMD analysis
sb0: DISTANCE ATOMS=1782,407  NOPBC # R217 - D129
sb1: DISTANCE ATOMS=1886,1300 NOPBC # R223 - D186
sb2: DISTANCE ATOMS=1948,531  NOPBC # R226 - D136
sb3: DISTANCE ATOMS=2004,407  NOPBC # R229 - D129
sb4: DISTANCE ATOMS=2057,1300 NOPBC # R232 - D186

MOVINGRESTRAINT ...
    ARG=rmsd1,rmsd2,sb0,sb1,sb2,sb3,sb4
    AT0={rmsd1s},{rmsd2s},{sb0s},{sb1s},{sb2s},{sb3s},{sb4s} STEP0=0         KAPPA0=0,0,0,0,0,0,0
    # turn on restraints over first 500 ps
    AT1={rmsd1s},{rmsd2s},{sb0s},{sb1s},{sb2s},{sb3s},{sb4s}  STEP1=125000    KAPPA1=500,500,500,500,500,500,500
    AT2={rmsd1e},{rmsd2e},{sb0e},{sb1e},{sb2e},{sb3e},{sb4e}   STEP2=500000    KAPPA2=500,500,500,500,500,500,500
    LABEL=smd
...

PRINT ARG=* FILE={output} STRIDE=2500"""

        with open(plumed, mode="w+") as f:
            f.writelines(SMD_STRING)


if __name__ == "__main__":
    ref = sys.argv[1:3]
    mode = sys.argv[3]
    plumed = sys.argv[4]
    output = sys.argv[5]
    if mode == "compute":
        cv_files = None
    elif mode == "pull":
        cv_files = sys.argv[6:8]
    else:
        raise ValueError("Unaccepted mode, must be either 'compute' or 'pull'")

    main(ref, mode, plumed, output, cv_files)
