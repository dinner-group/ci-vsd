set psfname civsd-pro
set dirxtc /project2/roux/scguo/ci-vsd/anton/xtc1000ns
set firstframe 0
set skipframe 1


set xtcname [lindex $argv 0]
set outfile /project2/roux/scguo/ci-vsd/data/projection_anton/$xtcname.txt

mol load psf ../$psfname.psf pdb ../$psfname.pdb
mol addfile $dirxtc/$xtcname.xtc waitfor all

set fil1 [open $outfile w]
puts $fil1

#cv delete
cv molid top
cv config "
    colvar {
    name distS4    
    distanceZ {
        main {
            psfSegID A
            atomNameResidueRange CA 217-233   # Select biased atoms from this file

            centerReference          # use relative coordinates
            rotateReference          # (translated and rotated frame of reference)
            refPositionsGroup {     
                atomsFile     civsd-pro-up.pdb   # Define frame of reference base on separate atom group from separate file
                atomsCol      B
                atomsColValue 2
            }
            refPositionsFile  civsd-pro-up.pdb   # initial coordinates for reference group
        }

        ref {
            dummyAtom (8.13, -3.96, 4.51)  # sets the zero position along the S4 axis; center of 217-233 CA
        }

        axis (-0.211, -0.282, 0.936)       # initial orientation is along the principle axis of S4 217-233 CA
   }
   }"
    
cv config "
    colvar {
    name spinS4
    spinAngle {
	componentCoeff    -1.0
        atoms {
            psfSegID A
            atomNameResidueRange CA 217-233

            centerReference          # use relative coordinates
            rotateReference          # (translated and rotated frame of reference)
            refPositionsGroup {      
                atomsFile     civsd-pro-up.pdb   # Define frame of reference base on separate atom group from separate file
                atomsCol      B
                atomsColValue 2
            }
            refPositionsFile  civsd-pro-up.pdb   # initial coordinates for reference group
        }
        refpositionsfile civsd-pro-up.pdb
    }
    }"

# animate delete all

for {set frame $firstframe} {$frame <= [molinfo top get numframes]} {incr frame $skipframe} {
    # animate read xtc $dirxtc/$xtcname.xtc beg $frame end $frame
    # animate goto $frame
    cv frame $frame
    cv update
    set cvs [cv printframe]
    set dist [lindex $cvs 1]
    set spin [lindex $cvs 2]
    puts $fil1 [format "%.3f %.3f" $dist $spin]
    # animate delete all
}
close $fil1
puts [molinfo top get numframes]
# mol delete all

exit

 
