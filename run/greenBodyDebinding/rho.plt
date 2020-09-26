#!/usr/bin/env gnuplot --persist -c

if (ARGC > 0) {
    set terminal postscript eps color font 'Helvetica,14'
    set output ARG1.".eps"
} else {
    set terminal qt size 1000, 600 font 'Helvetica,14'
}

set xlabel "Temperature (°C)"
set ylabel "Density of monomer (kg/m^3)"
set log y

KtoC(T) = T - 273.15

plot "postProcessing/minMax/0/fieldMinMax.dat" u (KtoC($4)):3 w l notitle lw 2

