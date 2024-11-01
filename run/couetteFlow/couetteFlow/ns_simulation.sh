#!/bin/bash

rm -rf ../ns
for U in $*; do
    mkdir -p ../ns/$U
    cp -r * ../ns/$U
    (
        echo "Simulate for U=$U"
        cd ../ns/$U
        rm -rf 0
        mv 0.ns 0
        U2=$(echo "print .5*$U" | python)
        sed -i~ "32s/1/$U2/" 0/U0
        ./Allrun
    )
done
