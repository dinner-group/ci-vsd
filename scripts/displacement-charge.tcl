# lappend auto_path {/home/shenr/bin/maplugin/tcllib1.5}
# package require math::statistics
package require pbctools

set psfname "/project/dinner/scguo/ci-vsd/models/MD-clustering-center/civsd.psf"
set dcdname [lindex $argv 0]
set outputname [lindex $argv 1]

set fil1 [open $outputname.dat w]

mol load psf $psfname
animate read dcd $dcdname

set firstframe 0
set lastframe [molinfo top get numframes]
set skipframe 1

for {set frame $firstframe} {$frame < $lastframe} {incr frame $skipframe} {
    # animate goto $frame
	set pbcbox [pbc get -first $frame -last $frame]
    puts $pbcbox
	set Lz [lindex [lindex $pbcbox 0] 2]
	#   protein
	set Qpro 0
	set sel1 [atomselect top "protein"]
	foreach i [$sel1 get z] j [$sel1 get charge] {
	    set Qpro [expr $Qpro + $j*($i+0.5*$Lz)/$Lz]
	}
	
	#   lipids
	set Qlip 0
	set sel2 [atomselect top "lipids"]
	foreach i [$sel2 get z] j [$sel2 get charge] {
	    set Qlip [expr $Qlip + $j*($i+0.5*$Lz)/$Lz]
	}
	
	#   water
	set Qwat 0
	set sel3 [atomselect top "water"]
	foreach i [$sel3 get z] j [$sel3 get charge] {
	    set Qwat [expr $Qwat + $j*($i+0.5*$Lz)/$Lz]
	}   
	
	#   ions
	set Qion 0
	set sel4 [atomselect top "ions"]
	foreach i [$sel4 get z] j [$sel4 get charge] {
	    if {$i>0} {
		set i [expr $i - $Lz]
	    }
	    set Qion [expr $Qion + $j*($i+1.0*$Lz)/$Lz] ;# unwrapped coordinates of ions 1.0*$Lz
	} 
	
	#   protein + lipids + water + ions	
	lappend Qd [expr $Qpro + $Qlip + $Qwat + $Qion]  

	$sel1 delete
	$sel2 delete
	$sel3 delete
	$sel4 delete	
	animate delete all
}
puts $fil1 $Qd
close $fil1

exit
