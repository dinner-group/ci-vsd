package require pbctools

set psfname civsd-all
set dcdname [lindex $argv 0]
set dirpsf /project/dinner/scguo/ci-vsd
set dirdcd /project/dinner/scguo/ci-vsd/pmepot


set firstframe [lindex $argv 1]
set lastframe [lindex $argv 2]

set molid [mol new $dirpsf/$psfname.prmtop]
mol addfile $dirdcd/$dcdname.dcd first $firstframe last $lastframe waitfor all

set iframe [molinfo top get frame]

# pbc unwrap -sel "protein"
pbc wrap -centersel "protein" -center com -compound residue -all

for {set i 0} {$i <= $iframe} {incr i} {
    set sel [atomselect top "protein and backbone" frame $i]
    foreach {xi yi zi} [measure center $sel] {}
    set pbcbox [pbc get -first $i -last $i]
    set px [lindex [lindex $pbcbox 0] 0]
    set py [lindex [lindex $pbcbox 0] 1]
    set pz [lindex [lindex $pbcbox 0] 2]
    set xm [expr $xi - 0.5*$px]
    set ym [expr $yi - 0.5*$py]
    set zm [expr $zi - 0.5*$pz]

    pbc wrap -sel "not protein" -compound res -shiftcenter "$xm $ym $zm" -first $i -last $i
    set sel1 [atomselect top all frame $i]
    set comx [expr -1 * $xi]
    set comy [expr -1 * $yi]
    set comz [expr -1 * $zi]
    $sel1 moveby "$comx $comy $comz"
}

animate write dcd $dirdcd/${dcdname}_wrap.dcd waitfor all

exit
