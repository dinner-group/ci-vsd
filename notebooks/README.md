# Notebooks

`figures` contains notebooks for loading data and rendering figures the paper.
  - See the `README.md` file in [figures](notebooks/figures) for more information
    
`anton2` contains the majority of calculations for the CVs, DGA quantities, and other quantities. It is organized as follows:
  - `cvs` contains notebooks for computing CVs (collective variables, e.g. salt-bridge distances) from the raw trajectories
  - `dga_calculations` contains notebooks for DGA calculation of kinetic statistics.
    The majority of statistics used for the figures are found in `dga_220325.ipynb`.
  - The remaining notebooks contain various other analyses, some of which are used for the figures.
    
`fitting` contains experiments for learning the committor using different methods (of which the LASSO was used in the paper).
