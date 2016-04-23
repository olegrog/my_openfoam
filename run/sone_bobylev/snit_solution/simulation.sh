#!/bin/bash

Kn=$(python -c "
import numpy, sys
#numpy.savetxt(sys.stdout, numpy.arange(0,0.07,0.0025), fmt='%.4f')
numpy.savetxt(sys.stdout, numpy.logspace(-2.6,-1,25), fmt='%.7f')
")

problem=$(basename $(pwd))
dir="../_$problem"

case=$1

[[ $# -eq 0 ]] && case=asym
echo "Use directory 0.$case"

rm -rf $dir
mkdir -p $dir
for kn in $Kn; do
    mkdir $dir/$kn
    cp -r * $dir/$kn/
    (
        echo "Simulate for Kn=$kn"
        cd $dir/$kn
        ./Allrun $case $kn
    )
done
