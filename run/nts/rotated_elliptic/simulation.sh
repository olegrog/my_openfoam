#!/bin/bash

for c in $*; do
    mkdir -p ../$c
    cp -r * ../$c/
    (
        echo "Simulate for $c"
        cd ../$c
        sed -i~ "4s/1/$c/" elliptic.geo
        ./Allrun
    )
done
