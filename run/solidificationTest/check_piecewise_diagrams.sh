#!/usr/bin/env bash

log="log.$(grep application system/controlDict | egrep -o '\w+Foam')"
[ -f "$log" ] || { echo "There is no $log file"; exit 1; }

gawk '
/^Selecting componentPhase/ {
    phase = $5
    name = $6
    phases[phase] = 0
}
/^ -- Ceq/ {
    T=$3; C=$6
    data[name][T][phase] = C
}
END {
    PROCINFO["sorted_in"] = "@ind_str_asc"
    for (name in data) {
        file = "_"name".txt"
        printf "T" > file
        for (phase in phases) {
            printf "\t%s\t", phase >> file
        }
        printf "\n" >> file
        for (T in data[name]) {
            printf "%.1f", T >> file
            for (phase in data[name][T]) {
                printf "\t%.3e", data[name][T][phase] >> file
            }
            printf "\n" >> file
        }
    }
}' $log

mapfile -t files < <(ls _*.txt)

gnuplot -p <<-EOF
    n = ${#files[@]}
    files = "${files[@]}"
    set key autotitle columnhead
    set multiplot layout n/2, 2
    set xlabel "Temperature (K)"
    #set format y "%.1e"
    do for [file in files] {
        stats file nooutput
        set title substr(file, 2, strstrt(file, '.') - 1)
        plot for [i=2:STATS_columns] file u 1:i w l
    }
EOF
