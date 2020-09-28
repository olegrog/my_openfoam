#!/usr/bin/env gnuplot --persist -c

file1 = "postProcessing/minMax/0/fieldMinMax.dat"
file2 = "postProcessing/minMaxKinetics/0/fieldMinMax.dat"
file3 = "data/strength.txt"

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
fromMPa(p) = p*1e6
filter(x,min,max) = (x > min && x < max) ? x : 1/0

unset label; set label "a)" at graph -0.11, 1
set ylabel "Mass fraction of the polymer component"
set yrange [0:1]
plot file2 u (KtoC($4)):6 w l title "Polymer1" lw 2, \
    file2 u (KtoC($4)):3 w l title "Polymer2" lw 2

unset label; set label "b)" at graph -0.15, 1
set ylabel "Diffusion coefficient (m^2/s)"
set log y
set yrange [1e-15:1e-5]
plot file1 u (KtoC($6)):8:9 w filledcurves notitle, \
    file1 u (KtoC($6)):8 w l lw 2 title "Min", \
    file1 u (KtoC($6)):9 w l lw 2 title "Max"

unset label; set label "c)" at graph -0.12, 1
set ylabel "Maximum density of monomer (kg/m^3)"
unset yrange
plot file1 u (KtoC($6)):3 w l notitle lw 2

unset label; set label "d)" at graph -0.15, 1
set ylabel "Pressure (Pa)"
set object 1 rect from  0,0 to 1,1
set object 1 rect fc rgb "cyan" fillstyle solid 1.0
plot file1 u (KtoC($6)):5 w l title "Maximum monomer pressure" lw 2, \
    file3 u 1:(fromMPa($2)) w lp title "Compressive strength", \
    file3 u 1:(fromMPa($3)) w lp title "Tensile strength", \
    file3 u 1:(fromMPa($4)) w lp title "Flexural strength", \
    file1 u (filter(KtoC($6), 100, 200)):5 with filledcurves x1 notitle fs transparent solid 0.25 lc rgb "black"

