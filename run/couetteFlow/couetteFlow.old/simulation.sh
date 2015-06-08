#!/bin/bash

for c in $*; do
    rm -rf ../$c
    mkdir -p ../$c
    cp -r * ../$c/
    (
        echo "Simulate for $c"
        cd ../$c
        wallU=$(echo "print $c*0.5" | python)
        sed -i~ "32s/1/$wallU/" 0/U0
        ./Allrun
    )
done
