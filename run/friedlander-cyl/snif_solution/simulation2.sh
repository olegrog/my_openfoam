#!/bin/bash

Kn=$(python -c "
import numpy, sys
numpy.savetxt(sys.stdout, numpy.logspace(-4,1,31), fmt='%.4f')
")

problem=$(basename $(pwd))
dir="../_$problem"

for kn in $Kn; do
    cp Allrun2 $dir/$kn/
    (
        echo "Simulate for Kn=$kn"
        cd $dir/$kn
        ./Allrun2 asym $kn
    )
done
