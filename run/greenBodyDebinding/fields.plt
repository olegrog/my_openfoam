#!/usr/bin/env -S gnuplot --persist -c

file1 = "postProcessing/minMax/0/fieldMinMax.dat"
file2 = "postProcessing/minMaxKinetics/0/fieldMinMax.dat"
file3 = "data/strength.txt"

if (ARGC > 0) {
    set terminal postscript eps color font 'Helvetica,11'
    set output ARG1.".eps"
    label_shift = 0.05
} else {
    set terminal qt size 1000, 700 font 'Helvetica,14'
    label_shift = 0
}

set encoding utf8
set multiplot layout 2, 2
set xlabel "Temperature (°C)"

KtoC(T) = T - 273.15
fromMPa(p) = p*1e6
filter(x,min,max) = (x > min && x < max) ? x : 1/0

unset label; set label "a)" at graph -0.11 - label_shift, 1
set ylabel "Mass fraction of a polymer component"
set yrange [0:1]
plot file2 u (KtoC($4)):6 w l title "Polymer1" lw 2, \
    file2 u (KtoC($4)):3 w l title "Polymer2" lw 2, \
    file1 u (KtoC($11)):13 w l title "Porosity" lw 2

unset label; set label "b)" at graph -0.15 - label_shift, 1
set ylabel "Mass diffusivity (m^2/s)"
set log y
set yrange [1e-15:1e-5]
set style fill transparent solid 0.5
plot file1 u (KtoC($11)):2:3 w filledcurves title "Due to convection (min/max)", \
    file1 u (KtoC($11)):6:7 w filledcurves title "Due to diffusion (min/max)", \
    file1 u (KtoC($11)):6 w l lw 2 lt -1 notitle, \
    file1 u (KtoC($11)):7 w l lw 2 lt -1 notitle, \
    file1 u (KtoC($11)):2 w l lw 2 lt -1 notitle, \
    file1 u (KtoC($11)):3 w l lw 2 lt -1 notitle, \

unset label; set label "c)" at graph -0.12 - label_shift, 1
set ylabel "Maximum density of monomer (kg/m^3)"
unset yrange
plot file1 u (KtoC($11)):5 w l notitle lw 2

unset label; set label "d)" at graph -0.15 - label_shift, 1
set ylabel "Pressure (Pa)"
set object 1 rect from  0,0 to 1,1
set object 1 rect fc rgb "cyan" fillstyle solid 1.0
set style fill transparent solid 1
plot file1 u (filter(KtoC($11), 100, 200)):9 with filledcurves x1 notitle lc rgb "gray",\
    file1 u (KtoC($11)):9 w l title "Maximum monomer pressure" lw 2, \
    file3 u 1:(fromMPa($2)) w lp title "Compressive strength", \
    file3 u 1:(fromMPa($3)) w lp title "Tensile strength", \
    file3 u 1:(fromMPa($4)) w lp title "Flexural strength", \

getValue(row,col,filename) = system('awk ''{if (NR == '.row.') print $'.col.'}'' '.filename.'')
stats file1 u 9 nooutput
print "Max pressure (MPa) = ", STATS_max/1e6, " when T (°C) = ", KtoC(getValue(STATS_index_max+3, 11, file1))

