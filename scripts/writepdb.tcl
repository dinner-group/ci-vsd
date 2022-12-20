mol load psf civsd-remd.psf namdbin clusterAAA/civsd.initial.BBB.coor
set sel [atomselect top all]
set x [expr XXX / 2.0]
set y [expr YYY / 2.0]
set z [expr ZZZ / 2.0]
$sel moveby "$x $y $z"
$sel writepdb clusterAAA/civsd.initial.BBB.pdb
exit
