#!/bin/bash

dir=postProcessing/singleGraph
[[ -d "$dir" ]] || { echo "Directory `$dir` is missing"; exit 1; }
dir="$dir/$(ls "$dir" | sort -g | tail -1)"
echo "Using $dir"

file0=constant/chemicalProperties
file1="$dir/line_lambda_p_rho.xy"
file2="$dir/line_U.xy"
file3=_solutionZND.txt

read_property() {
    local file=$1
    local prop=$2
    grep "^$prop " "$file" | sed 's/;.*//' | grep -Eo '[+-]?[0-9]+([.][0-9]+)?'
}

s=$(python3 <<PY
import numpy as np
x,u,v,w = np.loadtxt("$file2").T
idx = np.where(u**2>1e-6)
print(np.mean(x[idx[0][-2:]]))
PY
)
echo "Shock-wave coordinate = $s"

L=10
Q=$(read_property $file0 Q)
E=$(read_property $file0 E)
k=$(read_property $file0 k)

if [[ ! -f "$file3" ]]; then
    set -x
    ./find_ZND.py -v -L=$L -Q=$Q -E=$E -k=$k --output=$file3
    set +x
fi

gnuplot -p <<PLT
    set xrange [-$L:0]
    plot "$file1" u (\$1-$s):2 w l title "lambda", \
        "" u (\$1-$s):3 w l title "p", \
        "" u (\$1-$s):4 w l title "rho", \
        "" u (\$1-$s):(\$3/\$4) w l title "T", \
        "$file2" u (\$1-$s):2 w l title "U", \
        "$file3" u 1:2 w l lc 1 lw 2 dt 2 notitle, \
        "" u 1:3 w l lc 2 lw 2 dt 2 notitle, \
        "" u 1:4 w l lc 3 lw 2 dt 2 notitle, \
        "" u 1:(\$3/\$4) w l lc 4 lw 2 dt 2 notitle, \
        "" u 1:5 w l lc 5 lw 2 dt 2 notitle
PLT

