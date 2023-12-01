# Dynamics of activation in the voltage-sensing domain of Ci-VSP

This repository contains custom code for analyzing simulations of Ci-VSD for [ref. 1][1].

## Data
Example data to run the create [Figure 5](./notebooks/figures/FIG1_S4.ipynb) of the publication can
be obtained as a file `data.tar.gz` as Source Data with the accompanying publication. The tarball can be unzipped
using `tar -xzvf data.tar.gz`and the files used to generate the appropriate plots (see [notebooks](notebooks/figures/README.md)).

## Environment
Most of the analysis is performed in [Jupyter notebooks](./notebooks/). 
Notebooks to produce relevant figures are under the subfolder [figures](./notebooks/figures) along with a `README` file that provides further description. 
The data

The environment to setup the analysis is using Python 3.9.x, and detailed
packages used are listed in `requirements.txt`.
The analysis scripts have only been tested in 3.9.x and will not work in earlier versions (due to 
package dependencies).

### Dependencies
The primary dependencies are Jupyter notebooks
- Numerical analysis
    - `numpy<1.22`
    - `scipy`
    - `scikit-learn>=1.2`
    - `numba`
- Plotting
    - `matplotlib`
    - `seaborn`
    - `prettypyplot` ([https://braniii.gitlab.io/prettypyplot/](https://braniii.gitlab.io/prettypyplot/))
- MD analysis
    - `MDAnalysis>=2.0` ([MDAnalysis](https://www.mdanalysis.org/))
    - `pyemma` 2.5 ([PyEMMA](http://www.emma-project.org/))
    - `mdtraj` 1.9 ([MDTraj](https://www.mdtraj.org/))
- DGA/TPT analysis (from the Dinner group)
    - `extq` ([DGA](https://github.com/chatipat/extq))
    - `ivac` ([Integrated VAC](https://github.com/chatipat/ivac))

Other useful utility functions and plotting functions are found in 
`python/util.py` and `python/plotting.py`.

### DGA calculations
The code to run DGA calculations is provided in the `extq` package. The package can be installed
by cloning the repository `git clone https://github.com/chatipat/extq.git` and running `pip install -e extq`.

To reproduce the DGA calculations, run the notebook `dga.ipynb` in `notebooks/figures` using the [example
data](#data) provided.

### conda
The easiest way to set up is to use `requirements.txt` and create a new
Conda environment.
```
conda create --file requirements.txt -n civsd
```
Or one can create a new environment and install packages as necessary:
```
conda create -n civsd python=3.9
```

## Files
- Generic scripts to perform CV calculations are found in `scripts`.
These include several `tcl` scripts to use with [VMD](https://www.ks.uiuc.edu/Research/vmd/).
- Reference crystal structures (see [ref. 1][1]) and MD-simulation structures ([ref. 2][2]) are found in `models`.
- Scripts and starting files for setting up biased SMD (steered MD) simulations are in `biased/smd-prep`.
These were used to equilibrate new structures for unbiased simulations by pulling the protein forward and backwards.


## References
1. Guo, et al. [Dynamics of activation in the voltage-sensing domain of Ci-VSP.][1] *bioRxiv* **2022**.
1. Li, et al. [Structural mechanism of voltage-dependent gating in an isolated voltage-sensing domain.][2] *Nat. Struct. Mol. Bio.* **2014**.
1. Shen, et al. [Mechanism of Voltage Gating in the Voltage-Sensing Phosphatase Ci-VSP.][3]  *PNAS* **2022**. 

[1]: https://www.biorxiv.org/content/10.1101/2022.12.19.521128v2
[2]: http://www.nature.com/articles/nsmb.2768
[3]: https://www.pnas.org/doi/full/10.1073/pnas.2206649119
