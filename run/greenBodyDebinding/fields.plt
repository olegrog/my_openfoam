#!/usr/bin/env gnuplot --persist -c

file1 = "postProcessing/minMax/0/fieldMinMax.dat"
file2 = "postProcessing/minMaxKinetics/0/fieldMinMax.dat"

if (ARGC > 0) {
    set terminal postscript eps color font 'Helvetica,14'
    set output ARG1.".eps"
} else {
    set terminal qt size 1000, 700 font 'Helvetica,14'
}

set encoding utf8
set multiplot layout 2, 2
set xlabel "Temperature (°C)"

KtoC(T) = T - 273.15

set ylabel "Mass fraction of the polymer component"
set yrange [0:1]
plot file2 u (KtoC($4)):6 w l title "Polymer1" lw 2, \
    file2 u (KtoC($4)):3 w l title "Polymer2" lw 2

set ylabel "Diffusion coefficient (m^2/s)"
set log y
set yrange [1e-15:1e-5]
plot file1 u (KtoC($4)):6:7 w filledcurves notitle, \
    file1 u (KtoC($4)):6 w l lw 2 title "Min", \
    file1 u (KtoC($4)):7 w l lw 2 title "Max"

set ylabel "Maximum density of monomer (kg/m^3)"
unset yrange
plot file1 u (KtoC($4)):3 w l notitle lw 2

