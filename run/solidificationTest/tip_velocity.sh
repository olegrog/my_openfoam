#!/bin/bash

appl=$(grep application system/controlDict | egrep -o '\w+Foam')

gnuplot -p -e "set format y '%.1e'; set xrange [10:180]; plot '<&3' w l title 'Tip velocity'" \
    3< <(grep 'Tip v' log.$appl | awk '{print $4}')
