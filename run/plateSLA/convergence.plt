#!/usr/bin/env gnuplot --persist -c

if (ARGC > 0) {
    set terminal postscript eps color font 'Helvetica,14'
    set output ARG1.".eps"
} else {
    set terminal qt size 1000, 600 font 'Helvetica,14'
}

file = "log.polymerSolidFoam"
# Generate 3 columns
data = sprintf("<(awk '/for Dx/{a=$$8+0} /for Dy/{b=$$8+0} /for Dz/{c=$$8+0; print a,b,c}' %s)", file)

set xlabel "Iteration"
set ylabel "Residual"
set log y
set format y "%g"

plot data u 0:1 w l title 'Dx', '' u 0:2 w l title 'Dy', '' u 0:3 w l title 'Dz'
