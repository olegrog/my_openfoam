#!/bin/bash

./param.py

problem=$(basename $(pwd))
problem="../_$problem"

increase_ansemble() {
    pts=$(grep Part constant/dsmcProperties | sed 's/.*  *//;s/;.*//')
    pts=$(echo "print $pts/$1" | python)
    echo "print int(1/1.3806488e-23/$pts)" | python
    sed -i.orig "s/nEquivalentParticles.*/nEquivalentParticles            $pts;/" constant/dsmcProperties
}

for d in $(ls $problem | grep ^[0-9]); do
    [[ -f $problem/$d/log.nsfSimpleFoam ]] && continue
    (
        echo "Simulate for $d"
        cd $problem/$d
        #    increase_ansemble 4
        ./Allrun
    )
done
