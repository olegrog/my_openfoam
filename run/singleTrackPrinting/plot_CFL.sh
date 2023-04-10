#!/bin/bash

appl=$(grep application system/controlDict | egrep -o '\w+Foam')
maxCo=$(grep maxCo system/controlDict | awk '{print $2}' | sed s'/;//')
maxAlphaCo=$(grep maxAlphaCo system/controlDict | awk '{print $2}' | sed s'/;//')
maxDeltaTemp=$(grep maxDeltaTemp system/controlDict | awk '{print $2}' | sed s'/;//')

gnuplot -p -e "set yrange [0:2]; plot '<&3' w l title 'Courant',
    '<&4' w l title 'Interface', '<&5' w l title 'Temperature'" \
    3< <(awk '/^Courant Number/{print $6/'"$maxCo"'}' log.$appl) \
    4< <(awk '/^Interface Courant/{print $7/'"$maxAlphaCo"'}' log.$appl) \
    5< <(awk ' /^Temperature change/{print $6/'"$maxDeltaTemp"'}' log.$appl)
