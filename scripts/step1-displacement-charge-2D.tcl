#lappend auto_path {/u/sciteam/rshen/bin/maplugin/tcllib1.5}

lappend auto_path {/home/shenr/bin/maplugin/tcllib1.5}
package require math::statistics
package require pbctools

set psfname civsd-remd
set dirpsf ../
set cvsfile change_wins_replicas
set outputname civsd-gating-charges-neutral-2d

set wins [exec awk {{printf("%.0f ", $1)}} ./$cvsfile.txt ];
set dists [exec awk {{printf("%.2f ", $3)}} ./$cvsfile.txt ];
set rots [exec awk {{printf("%.1f ", $4)}} ./$cvsfile.txt ];

set fil1 [open ./data/$outputname.dat w]

mol load psf $dirpsf/$psfname.psf
foreach win $wins dist $dists rot $rots {  
    set dirdcd ../output/$win
    set firstframe 100
    set lastframe 600
    set skipframe  10
    set dcdname civsd.job4.$win.sort
    set Qd {}
    
    for {set frame $firstframe} {$frame < $lastframe} {incr frame $skipframe} {
	animate read dcd $dirdcd/$dcdname.dcd beg $frame end $frame
	set pbcbox [pbc get -now]
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

    set Q_mean [::math::statistics::mean $Qd]
    set Q_std [::math::statistics::stdev $Qd]
    
    set output_dis [format "%.0f %.2f %.1f %.3f %.3f" $win $dist $rot $Q_mean $Q_std]
    puts $fil1 $output_dis
}
close $fil1
mol delete all

exit

