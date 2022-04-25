#!/bin/bash

log=log.detonationFoam
[ -f "$log" ] || { echo "There is no $log file"; exit 1; }

dict=system/movingFrameDict
[ -f "$dict" ] || { echo "There is no $dict file"; exit 1; }

if [ -z $DISPLAY ]; then
    cmd='set term dumb;'
fi

gnuplot -p -e "$cmd plot '<&3' u 2:1 w l title 'Urel(t)'" \
    3< <(awk '/^Time =/{print $3} /^Urel =/{printf $3" "}' $log | sed 's/(//')

if grep -q 'adaptive *1;' system/controlDict; then
    gnuplot -p -e "$cmd plot '<&3' u 2:1 w l title '<lambda>(t)'" \
    3< <(awk '/^Time =/{print $3} /^meanValue =/{printf $3" "}' $log)
fi

