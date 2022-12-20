set pdb [lindex $argv 0]
mol addfile $pdb

for {set i 0} {$i < 24} {incr i} {
    [atomselect top all frame $i] writepdb $i.pdb
}

exit
