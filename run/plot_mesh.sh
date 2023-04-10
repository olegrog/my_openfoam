#!/bin/bash -e

if [ -f system/controlDict ]; then
    log="log.$(grep application system/controlDict | egrep -o '\w+Foam')"
    [ -f "$log" ] || { echo "There is no $log file"; exit 1; }
else
    log="$(ls log.*Foam)"
fi

if [ -z $DISPLAY ]; then
    cmd='set term dumb;'
fi

gnuplot -p -e "$cmd; set title 'Number of cells'; plot '<&3' w l notitle" \
    3< <(grep 'cells\.' $log | awk '{print $5}')
