import sys
import os
import mdtraj as md
from mdtraj import Trajectory
from mdtraj.formats import DTRTrajectoryFile
from mdtraj.utils import in_units_of

sys.path.insert(1, "/project/dinner/scguo/ci-vsd/python")
from util import zero_h_occupancy

TOPFILE = "/project/dinner/scguo/ci-vsd/unbiased/000/civsd.prmtop"


def save_frame(trajfile, frame, outfile):
    if not trajfile.endswith("dtr"):
        fr = md.load(
            trajfile,
            top=TOPFILE,
            frame=frame,
        )
    elif trajfile.endswith("dtr"):
        # for Anton trajectories
        traj = DTRTrajectoryFile(trajfile)
        traj.seek(frame)
        xyz, time, box_length, box_angle = traj.read(n_frames=1)
        # default Trajectory units are nm, need to convert from Angstroms
        in_units_of(xyz, 'angstroms', Trajectory._distance_unit, inplace=True)
        in_units_of(box_length, 'angstroms', Trajectory._distance_unit, inplace=True)

        # create new trajectory
        top = md.load_topology(TOPFILE)
        fr = Trajectory(xyz=xyz, topology=top, time=time,
                        unitcell_lengths=box_length,
                        unitcell_angles=box_angle)
    fr.save_pdb(outfile)


def strip_water(pdbfile, outfile=None):
    if outfile is None:
        outfile = "protein.pdb"
    struct = md.load(
        pdbfile,
        top=TOPFILE,
    )
    struct_prot = struct.atom_slice(struct.topology.select("protein"))
    struct_prot.save_pdb(outfile)
    return outfile


def main(trajfile, frame, outfolder):
    outfile = f"{outfolder}/civsd.pdb"
    save_frame(trajfile, frame, outfile)
    nowat = f"{outfolder}/civsd-nowat.pdb"
    _ = strip_water(outfile, nowat)
    zero_h_occupancy(nowat, new_file_name=nowat, savelines=3)


if __name__ == "__main__":
    trajfile = sys.argv[1]
    frame = int(sys.argv[2])
    outfolder = sys.argv[3]
    if not os.path.exists(outfolder):
        os.mkdir(outfolder)

    main(trajfile, frame, outfolder)
