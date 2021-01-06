#!/bin/bash

appl=$(grep application system/controlDict | egrep -o '\w+Foam')

gnuplot -p -e "set format y '%.1e'; plot '<&3' w l" \
    3< <(grep absorp log.$appl | awk '{print $13}')
