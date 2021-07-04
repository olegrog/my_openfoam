#!/bin/bash -e

if [ -f system/controlDict ]; then
    log="log.$(grep application system/controlDict | egrep -o '\w+Foam')"
    [ -f "$log" ] || { echo "There is no $log file"; exit 1; }
else
    log="$(ls log.*Foam)"
fi


gnuplot -p -e "set format y '%.1e'; set log y; plot '<&3' w lp" \
    3< <(grep deltaT $log | awk '{print $3}')
