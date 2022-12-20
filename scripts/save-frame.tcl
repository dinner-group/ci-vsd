# Saves a single frame from a trajectory as a pdb file
# Takes the trajectory file (xtc) as the first argument, the desired frame number
# as the second argument, and the name of the output file as the 3rd.

# load trajectory
mol new /project2/roux/scguo/ci-vsd/civsd-pro.pdb
set traj [lindex $argv 0]
set desiredframe [lindex $argv 1]
mol addfile "${traj}" waitfor all
set out [lindex $argv 2]

# Write desired frame to pdb file
set selection [atomselect top "all"]
$selection frame $desiredframe 
$selection update
$selection writepdb $out

exit
