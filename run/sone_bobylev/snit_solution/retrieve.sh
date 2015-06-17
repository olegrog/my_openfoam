#!/bin/bash

. $WM_PROJECT_DIR/bin/tools/RunFunctions

problem=$(basename $(pwd))
dir="../_$problem"

cd $dir

printf '#%11s %12s\n' Kn "int(T_y=0.5)"
for Kn in $(ls | sort -g); do
    intT=$(grep 'over area' $Kn/log.patchIntegrate | tail -1 | awk '{ print $11 }')
    printf "%.6e %.5e\n" $Kn $intT
done
