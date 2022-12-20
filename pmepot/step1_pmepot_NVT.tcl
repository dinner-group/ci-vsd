
package require pmepot

set psfname civsd-all
set dcdname [lindex $argv 0]
set dirpsf /project/dinner/scguo/ci-vsd
set dirdcd /project/dinner/scguo/ci-vsd/pmepot


set firstframe [lindex $argv 1]
set lastframe [lindex $argv 2]

set molid [mol new $dirpsf/$psfname.prmtop]
mol addfile $dirdcd/$dcdname.dcd first $firstframe last $lastframe waitfor all

set iframe [molinfo $molid get frame]
set sel [atomselect $molid "protein"]

# cell: {0 0 0} {76.7581525578 0 0} {0 76.7581525578 0} {0 0 91.4578865745}
#pmepot -mol $molid -frames all -ewaldfactor 0.25 -cell {{0 0 0} {103 0 0} {0 103 0} {0 0 106}} -grid {220 220 240} -dxfile prod_dmpcps.dx 
#pmepot -mol $molid -frames all -ewaldfactor 0.25 -grid {100 100 100} -dxfile prod_popc.dx -xscfile ./dcdfile/civsd-fep-1.0.xsc

pmepot -mol $molid -frames all -ewaldfactor 0.25 -grid 0.5 -dxfile $dcdname.dx 
# pmepot -sel $sel -frames all -ewaldfactor 0.25 -grid 0.5 -dxfile $dcdname$firstframe.dx 


exit
