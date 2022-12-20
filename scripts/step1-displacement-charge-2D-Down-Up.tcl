#lappend auto_path {/u/sciteam/rshen/bin/maplugin/tcllib1.5}

# Down: win68  (-4.0 -45.0)
# Up:   win102 ( 0.5  -5.0)
# residues: 106-239

lappend auto_path {/home/shenr/bin/maplugin/tcllib1.5}
package require math::statistics
package require pbctools

set psfname civsd-remd
set dirpsf ../
set outputname civsd-gating-charges-neutral-2d-down-up

if { [info exists elist] } {
	unset elist
} 

if { [info exists Qd] } {
	unset Qd
}  
	
	
set fil [open ./data/$outputname.txt w]

foreach win {68 102} {

	mol load psf $dirpsf/$psfname.psf
    set dirdcd ../output/$win
    set firstframe 100
    set lastframe 600
    set skipframe  10
    set dcdname civsd.job4.$win.sort
      
    	
	for {set frame $firstframe} {$frame < $lastframe} {incr frame $skipframe} {
		animate read dcd $dirdcd/$dcdname.dcd beg $frame end $frame
		set pbcbox [pbc get -now]
		set Lz [lindex [lindex $pbcbox 0] 2]
		#   protein
		
		for {set i 106} {$i<=239} {incr i} {
		
			set Q1 0
			set sel1 [atomselect top "protein and resid $i"]
			foreach ii [$sel1 get z] jj [$sel1 get charge] {
	    		set Q1 [expr $Q1 + $jj*($ii+0.5*$Lz)/$Lz]
			}
			lappend Qd($win,$i) $Q1
			$sel1 delete
		}	
		animate delete all
	}


	for {set i 106} {$i<=239} {incr i} {   
	
    	set Q_mean [::math::statistics::mean $Qd($win,$i)]
    	#set Q_std [::math::statistics::stdev $Qd($win,$i)]
    	
    	lappend elist($win,$i) $Q_mean
    }
}

set delQ 0

for {set i 106} {$i<=239} {incr i} {
    set delQ [expr $delQ + ($elist(102,$i)-$elist(68,$i))]
    set output_dis [format "%.0f %.3f %.3f %.3f %.3f" $i $elist(68,$i) $elist(102,$i) [expr $elist(102,$i)-$elist(68,$i)] $delQ]
    puts $fil $output_dis
}

close $fil

mol delete all

exit

