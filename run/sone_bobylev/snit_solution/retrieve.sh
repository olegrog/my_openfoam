#!/bin/bash

. $WM_PROJECT_DIR/bin/tools/RunFunctions

problem=$(basename $(pwd))
dir="../_$problem$1"
kn_corr=$HOME/latex/science/7_sone_bobylev/contours/knudsen_layer_correction.py

cd $dir

get_value() {
    macro=$1
    patch=$2
    grep 'over area' $Kn/log.patchIntegrate.$macro.$patch | tail -1 | awk '{ print $11 }' | egrep -o '[0-9e\.\-]+' | awk '{ print $1*10 }'
}

printf '#%11s %11s %12s %11s %11s\n' Kn topT topU bottomT bottomU
for Kn in $(ls | sort -g); do
    last=$(ls $Kn | grep '^[0-9]' | tail -1)
    $kn_corr $Kn/VTK/$Kn'_'$last.vtk $Kn/tmp.vtk $Kn
done
