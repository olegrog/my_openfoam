#!/bin/bash

problem=$(basename $(pwd))
cd "../_$problem$1"

get_value() {
    macro=$1
    patch=$2
    grep 'over area' $dir/log.patchIntegrate.$macro.$patch | tail -1 | awk '{ print $11 }' | egrep -o '[0-9e\.\-]+' | awk '{ print $1 }'
}

printf '#%11s %11s %11s\n' Kn delta_p2 p0
for dir in $(ls | sort -g); do
    Kn=$(head -1 $dir/*.geo | awk '{ print $3 }' | sed 's/;//')
    [[ -f $dir/log.calcGradients ]] || continue
    last=$(grep Time $dir/log.calcGradients | tail -1 | awk '{ print $3 }')
    p0=$(grep 'internalField' $dir/$last/p0 | awk '{ print $3 }' | sed 's/;//')
    delta_p2=$(echo "print '%.5e' % (4*(($(get_value p2 right))-($(get_value p2 left)))/$p0)" | python)
    printf "%.6e %.5e %.5e\n" $Kn $delta_p2 $p0
done
