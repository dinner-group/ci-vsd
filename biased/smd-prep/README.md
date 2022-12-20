# Directory
run01 : down-seed2, testing run
run02 : down-seed2, testing run 2
run03 : down-seed2, longer run (2ns) w/ more on z-axis
run04 : down-seed2, 2 ns run with projection CV from abmd/run_10/
run05 : down-seed2, 2 ns run with RMSD to up-seed2.pdb, plumed-2 with higher force constant (2000)
run06 : down-seed2, 2 ns run with DRMSD to up-seed2.pdb
run07 : down-seed2, 4 ns run with RMSD to upplus-seed.pdb, k = 2000
run08 : down-seed2, 4 ns run with 6 contacts (rational definition)
run09 : down-seed2, 10 ns run with RMSD to down and upplus, 6 salt bridges from ABMD (Cz)
run10 : down-seed2, 10 ns run with RMSD to down and upplus, 6 salt bridges from ABMD (Cz), k = 3000
run11 : run10 rst, 10 ns run with RMSD from start and down-seed2 states, 6 salt bridges from ABMD (Cz)
run12 : run11 rst, 10 ns run with RMSD from start and up-seed2 states, 6 salt bridges from ABMD (Cz)
run13 : run12 rst, 10 ns run with RMSD from start and down-seed2 states, 6 salt bridges, k = 200 (Rong's settings)
run14 : down-seed2, 10 ns run with RMSD from start and up-seed2 states, k = 200 (Rong's settings)
run15 : down-seed2, 10 ns run with RMSD from start and downminus-seed.pdb, k = 200
run16 : up-seed2, 10 ns run with RMSD from start and upplus-seed.pdb, k = 200
run17 : down-seed2, 10 ns run with heavy atom RMSD from start and up-seed2.pdb, k = 200/500
run18 : down-seed2, 10 ns run with heavy atom RMSD from start and downminus-seed.pdb, k = 200/500
run19 : up-seed2, 10 ns run with heavy atom RMSD from start and upplus-seed.pdb, k = 200/500/1000
run20 : run13 rst, 10 ns run with RMSD from start and up-seed2 states, 6 salt bridges, k = 200
run21 : run21 rst, 10 ns run with RMSD from start and down-seed2 states, 6 salt bridges, k = 200
run22 : run21 rst, 10 ns run with RMSD from start, 6 salt bridges, ../../unbiased/085/civsd-4.nc 27, k = 500
run23 : run21 rst, 10 ns run with RMSD from start, 6 salt bridges, ./../unbiased/091/civsd-2.nc 2429, k = 500
run24 : run21 rst, 10 ns run with RMSD from start, 6 salt bridges, ./../unbiased/293/civsd-2.nc 1084, k = 500
run25 : run21 rst, 10 ns run with RMSD from start, 6 salt bridges, ./../unbiased/115/civsd.nc 1103, k = 500
run26 : run21 rst, 10 ns run with RMSD from start, 6 salt bridges, ./../unbiased/102/civsd-2.nc 52, k = 500
run27 : run21 rst, 10 ns run with RMSD from start, 6 salt bridges, ./../unbiased/284/civsd.nc 2907, k = 500

