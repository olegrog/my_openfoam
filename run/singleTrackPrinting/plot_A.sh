#!/bin/bash

print_help() {
    cat << EOF
Plot total absorptivity versus time.

Usage: ./$(basename "$0") [<options>]
Options:
  -min=<value>          Set a minimum A value
  -max=<value>          Set a maximum A value
  -eps                  Create a EPS file instead of graphical terminal
  -help                 Print this help
EOF
    exit 1
}

appl=$(grep application system/controlDict | egrep -o '\w+Foam')
file=A.eps
Amin=0
Amax=1

for arg; do case $arg in
    -min=*)     Amin="${arg#*=}";;
    -max=*)     Amax="${arg#*=}";;
    -eps)       eps=1;;
    -h|-help)   print_help;;
    -*)         echo "Unknown option '$arg'."; print_help;;
    *)          echo "Unknown argument '$arg'."; print_help;;
esac; done


if [[ -n "$eps" ]]; then
    cmd="set terminal postscript eps color font 'Helvetica,14'; set output '$file';"
fi

gnuplot -p -e "$cmd set yrange [$Amin:$Amax]; plot '<&3' w l notitle" \
    3< <(awk '/^Time =/{printf $3" "}/effective absorptivity/{print $13}' log.$appl)

if [ -f "$file" ] && command -v epstopdf > /dev/null; then
    epstopdf $file
fi
