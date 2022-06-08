#!/bin/bash

print_help() {
    cat << EOF
Plot laser properties versus time.

Usage: ./$(basename "$0") [<options>]
Options:
  -eps          Create a EPS file instead of graphical terminal
  -quiet        Do not print filenames that are being processed
  -log          Read a log file instead
  -help         Print this help
EOF
    exit 1
}

appl=$(grep application system/controlDict | egrep -o '\w+Foam')
output=laser.eps
tmpfile=_laser.txt

for arg; do case $arg in
    -eps)           eps=1;;
    -h|-help)       print_help;;
    -q|-quiet)      quiet=1;;
    -l|-log)        log=log.$appl;;
    -l=*|--log=*)   log="${arg#*=}";;
    -*)             echo "Unknown option '$arg'."; print_help;;
    *)              echo "Unknown argument '$arg'."; print_help;;
esac; done


if [[ -n "$eps" ]]; then
    cmd="set terminal postscript eps color font 'Helvetica,14'; set output '$output';"
fi

get_scalar() {
    local file=$1
    local property=$2
    grep "$property" "$file" | awk '{print substr($2, 1, length($2)-1)}'
}

get_vector1() {
    local file=$1
    local property=$2
    grep "$property" "$file" | awk '{print $3}'
}

[[ -f "$tmpfile" ]] && echo "File $tmpfile already exists."
if [[ "$log" ]]; then
    [[ -f "$log" ]] || { echo "File $log does not exist."; exit 1; }
    [[ $quiet ]] || echo "Processing file $log..."
fi

if [[ $log ]]; then
    awk '
        /^Time =/ { printf $3" " }
        /^Laser: / { print $4,$7,substr($10,2),substr($15,2) }
    ' "$log" > $tmpfile
else
    for f in $(find . -name uniform | sed s_\./__ | grep ^[0-9] | sort -g); do
        time=$(get_scalar $f/time value)
        [[ -f $tmpfile ]] && grep -Fq $time $tmpfile && continue
        power=$(get_scalar $f/laser power)
        radius=$(get_scalar $f/laser radius)
        velocity=$(get_vector1 $f/laser velocity)
        position=$(get_vector1 $f/laser position)
        switchedOn=$(get_scalar $f/laser switchedOn)
        power=$(bc <<< "$power*$switchedOn")
        echo "$time $power $radius $velocity $position" >> $tmpfile
        [[ $quiet ]] || echo $f
    done
fi

gnuplot -p <<PLT
    $cmd
    set multiplot layout 2,2
    plot "$tmpfile" u 1:2 w l title "power"
    plot "$tmpfile" u 1:3 w l title "radius"
    plot "$tmpfile" u 1:4 w l title "velocity"
    plot "$tmpfile" u 1:5 w l title "position"
PLT

if [ -f "$output" ] && command -v epstopdf > /dev/null; then
    epstopdf $output
fi
