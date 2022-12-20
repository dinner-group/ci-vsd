import os
import sys
import mdtraj as md

sys.path.insert(1, "/project/dinner/scguo/ci-vsd/python")
from util import zero_h_occupancy
from prep import strip_water


def main(input, output):
    _ = strip_water(input, output)
    zero_h_occupancy(output, new_file_name=output, savelines=3)


if __name__ == "__main__":
    input, output = sys.argv[1], sys.argv[2]
    main(input, output)
