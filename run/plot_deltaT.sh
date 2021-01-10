#!/bin/bash

appl=$(grep application system/controlDict | egrep -o '\w+Foam')

gnuplot -p -e "set format y '%.1e'; set log y; plot '<&3' w lp" \
    3< <(grep deltaT log.$appl | awk '{print $3}')
