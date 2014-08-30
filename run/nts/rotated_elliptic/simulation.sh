#!/bin/bash

for c in $*; do
    mkdir -p ../$c
    cp -r * ../$c/
    (
        echo "Simulate for $c"
        cd ../$c
        sed -i~ "5s/1/$c/" elliptic-transfinite.geo
        ./Allrun
    )
done
