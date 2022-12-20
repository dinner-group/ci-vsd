package require hbonds

# trajectory file is 1st argument
set traj [lindex $argv 0]
# output file is 2nd argument
set out [lindex $argv 1]

# load
mol new "/project2/roux/scguo/ci-vsd/unbiased/000/civsd.pdb"
mol addfile $traj waitfor all

# find HBonds between R217 and phosphate head groups
set phos [atomselect top  \
    "name O11 O12 O13 O14 and lipid and within 4 of resid 217"]
set arg [atomselect top "resid 217"]

# update selection each frame
hbonds -sel1 $phos \
    -sel2 $arg \
    -dist 4.0 \
    -ang 30 \
    -upsel "yes" \
    -plot no \
    -writefile yes \
    -outfile $out \
    -type all
