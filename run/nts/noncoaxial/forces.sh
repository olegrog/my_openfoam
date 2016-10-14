#!/bin/bash

printf "%7s %11s %11s %11s %11s %11s %11s\n" '# alpha' p0 inner outer hydro thermal visco
for d in $(ls | grep '[0-9]' | sort -n); do
    tau=$(echo "$d-1" | bc)
    last=$(grep '^Time =' $d/log.snitSimpleFoam | tail -1 | awk '{ print $3 }')
    p0=$(grep uniform $d/$last/p0 | awk '{ print $3 }' | sed 's/;//')
    inner=$(tail $d/log.wallSnitForces | grep inner | awk '{ print $2 }' | sed 's/(//')
    outer=$(tail $d/log.wallSnitForces | grep outer | awk '{ print $2 }' | sed 's/(//')
    hydro=$(tail $d/log.wallSnitForces | grep hydro | awk '{ print $2 }' | xargs | sed 's/(//g')
    visco=$(tail $d/log.wallSnitForces | grep hydro | awk '{ print $6 }' | xargs | sed 's/(//g')
    thermal=$(tail $d/log.wallSnitForces | grep hydro | awk '{ print $10 }' | xargs | sed 's/(//g')
    printf "%.1e % .4e % .4e % .4e % .4e % .4e % .4e % .4e % .4e % .4e\n" $tau $p0 $inner $outer $hydro $thermal $visco
done
