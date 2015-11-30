#!/bin/bash

get_p2() {
    grep 'area magnitude' log.patchIntegrate.p2.$1 | awk '{ print $11 }'
}

get_p0() {
    for t in $(grep 'Time =' log.patchIntegrate.p2.left | awk '{ print $3 }'); do
        grep 'internalField' $t/p0 | awk '{ print $3 }' | sed 's/;//'
    done
}

left=$(mktemp)
right=$(mktemp)
press=$(mktemp)

get_p2 left > $left
get_p2 right > $right
get_p0 > $press

paste $left $right $press | awk '{ print 100*4*($2-$1)/$3 }'

rm $left $right $press
