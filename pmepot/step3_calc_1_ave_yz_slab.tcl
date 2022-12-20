# RSHEN; 5-24-2019; 10-18-2020
#
# neg => delV (-); pos => delV (+)
# 
package require math::statistics


set sysname "R223A"
set confs "down-minus"

# offsets from "calc_0_integrate_z_axis.tcl"
#

set offsets "-0.1826 -0.2056" ;# R223A

#
#  0.5 V for neg; 
# -0.5 V for pos
#

set mempots "negative positive"
set delVs "0.5 -0.5"

foreach mempot $mempots delV $delVs offset $offsets {

	set inpname civsd-$confs-gc-$sysname-$mempot-500mV-xyz2.5
	set outname civsd-$confs-gc-$sysname-$mempot-500mV-yz-bias

	if { [info exists elist] } {
		unset elist
	} 

	set y {}
	set z {}

	set infile [open ./data/$inpname.txt r]
	
	while {[gets $infile line]>=0} {
		if {[llength $line]>0} {
			set y0 [lindex $line 0]
			lappend y $y0
			set z0 [lindex $line 1]
			lappend z $z0
			set ene [lindex $line 2]
			lappend elist($y0,$z0) $ene
		}
	}
	close $infile

	set outfile [open ./data/$outname.txt w]
	set ylist [lsort -unique -real -increasing $y]
	set zlist [lsort -unique -real -increasing $z]

	set num [expr [llength $zlist] -1]
	set del [expr $delV / $num]

	foreach yy $ylist {
  	  set j 0
		foreach zz $zlist {
			set ene_mean  [expr [::math::statistics::mean $elist($yy,$zz)] - $offset]
			set ene_mean_corr [expr $ene_mean + $j*$del]
			
			if {$ene_mean_corr > -0.05 && $ene_mean_corr < 0.0} {
				set ene_mean_corr 0.0
			}
			
			puts $outfile [format "%12.3f %12.3f %12.3f %12.3f" $yy $zz $ene_mean $ene_mean_corr]
			incr j
		}
	}
	close $outfile
}

exit