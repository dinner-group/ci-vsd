mol new /project/dinner/scguo/ci-vsd/civsd-pro.psf
mol addfile [lindex argv 0]
package require namdenergy

set r226 [atomselect top "protein and resid 226"]
set d129 [atomselect top "protein and resid 129"]
namdenergy -vdw -elec -sel $r226 $d129 -ofile "sb0.csv" -extsys step7_production.restart.xsc -exe "namd2_multicore" -par toppar/par_all36_carb.prm -par toppar/par_all36_cgenff.prm -par toppar/par_all36_lipid.prm -par toppar/par_all36_na.prm -par toppar/par_all36m_prot.prm -par toppar/toppar_water_ions.str -par ../unk/unk.prm
