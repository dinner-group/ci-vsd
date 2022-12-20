# On midway

module load amber/16
module load vmd

wins=0
for (( i=0; i < 180; i++ ))
do
    del=$(awk "BEGIN {printf \"%.0f\n\", $i%3}")
    if [ $del = 0 ]
    then
        mkdir cluster$wins
        cp ./input/civsd.initial.$i.coor ./cluster$wins
        cp ./input/civsd.initial.$i.xsc ./cluster$wins
        xxx=$(awk '($1!="#" && $1!="#$LABELS") {printf("%.1f ", $2)}' ./cluster$wins/civsd.initial.$i.xsc)
        yyy=$(awk '($1!="#" && $1!="#$LABELS") {printf("%.1f ", $6)}' ./cluster$wins/civsd.initial.$i.xsc)
        zzz=$(awk '($1!="#" && $1!="#$LABELS") {printf("%.1f ", $10)}' ./cluster$wins/civsd.initial.$i.xsc)
        sed -e "s/AAA/$wins/g" -e "s/BBB/$i/g"  -e "s/XXX/$xxx/g" -e "s/YYY/$yyy/g" -e "s/ZZZ/$zzz/g" < writepdb.tcl > writepdb-x.tcl
        vmd -dispdev text -e writepdb-x.tcl

        echo "$xxx $yyy $zzz"
        chamber -xpsf civsd-remd.psf -crd ./cluster$wins/civsd.initial.$i.pdb -top topparC36/top_all36_prot.rtf -param topparC36/par_all36_prot.prm -str topparC36/top_all36_lipid.rtf topparC36/top_water_ions.inp topparC36/par_all36_lipid.prm topparC36/par_water_ions_nbfix.prm -p ./cluster$wins/civsd_$wins.prmtop -inpcrd ./cluster$wins/civsd_$wins.inpcrd -cmap -box $xxx $yyy $zzz -vmd

        sed -e "s/AAA/$wins/g" < amber-hmr.parmed > amber-hmr-x.parmed
        parmed.py ./cluster$wins/civsd_$wins.prmtop amber-hmr-x.parmed

        let wins="$wins+1"
    fi
done



