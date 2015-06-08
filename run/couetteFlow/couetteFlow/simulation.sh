#!/bin/bash

Kn='
3.16228e-03
3.98107e-03
5.01187e-03
6.30957e-03
7.94328e-03
1.00000e-02
1.25893e-02
1.58489e-02
1.99526e-02
2.51189e-02
3.16228e-02
3.98107e-02
5.01187e-02
6.30957e-02
7.94328e-02
1.00000e-01
1.25893e-01
1.58489e-01
1.99526e-01
2.51189e-01
3.16228e-01
'

for U in $*; do
    rm -rf ../$U
    mkdir -p ../$U
    for kn in $Kn; do
        mkdir -p ../$U/$kn
        cp -r * ../$U/$kn/
        (
            echo "Simulate for U=$U, Kn=$kn"
            cd ../$U/$kn
            sed -i~ "37s/2/$U/" 0/U0
            sed -i~ "28s/1/$kn/" constant/transportProperties
            ./Allrun
        )
    done
done
