#!/bin/bash

problem=$(basename $(pwd))
cd "../_$problem$1"

printf '#%11s %11s %11s %11s %11s\n' Kn maxU p0 maxUx maxUy
for dir in $(ls | sort -g); do
    Kn=$(grep kn $dir/constant/nondimensionalProperties | awk '{ print $2 }' | sed 's/;//')
    [[ -d $dir/postProcessing ]] || continue
    maxU=$(tail -1 $dir/postProcessing/interiorMax/0/cellSource.dat | awk '{ print $3 }')
    #maxU=$(tail -1 $dir/postProcessing/minMax/0/fieldMinMax.dat | awk '{ print $7 }')
    maxUx=$(tail -1 $dir/postProcessing/minMax/0/fieldMinMax.dat | awk '{ print $8 }' | sed 's/(//')
    maxUy=$(tail -1 $dir/postProcessing/minMax/0/fieldMinMax.dat | awk '{ print $9 }')
    last=$(grep '^Time =' $dir/log.snitSimpleFoam | tail -1 | awk '{ print $3 }')
    p0=$(grep 'internalField' $dir/$last/p0 | awk '{ print $3 }' | sed 's/;//')
    printf "%.6e %.5e %.5e %.5e %.5e\n" $Kn $maxU $p0 $maxUx $maxUy
done
