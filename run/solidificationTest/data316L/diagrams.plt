#!/usr/bin/env gnuplot -persist
EPS=0

if (EPS) {
    set terminal postscript eps color font 'Helvetica,10' size 5in, 3in
    set output "1.eps"
} else {
    set term qt 0 position 200, 200
}

set multiplot layout 2, 2
set key autotitle columnhead
set xlabel "Temperature (K)"

DT = 273.15
TS = 1698.25
TS2 = 1707.15
TA = 1714.23
TL = 1726.78
C = 0.03; Cr = 0.17; Mo = 2.5; Ni = 0.12

rhoF_(T, C, Cr, Mo, Ni) = 7875.96 - 0.2970*T - 5.62e-5*T**2 \
    + (-206.35 + 0.00778*T + 1.472e-6*T**2)*C \
    + (-8.58 + 1.229e-3*T + 0.852e-7*T**2 + 0.018367*Cr)*Cr \
    + 30.78*Mo \
    + (-0.22 - 0.470e-3*T - 1.855e-7*T**2 + 0.104608*Ni)*Ni
rhoA_(T, C, Cr, Mo, Ni) = 8099.79 - 0.5060*T \
    + (-118.26 + 0.00739*T)*C \
    + (-7.59 + 3.422e-3*T - 5.388e-7*T**2 - 0.014271*Cr)*Cr \
    + 12.45*Mo \
    + (1.54 + 2.267e-3*T - 11.26e-7*T**2 + 0.062642*Ni)*Ni
rhoL2(T, C, Cr, Mo, Ni) = 8319.49 - 0.835*T \
    + (-83.19 + 0.00835*T)*C \
    + (-14.77 + 0.00535*T)*Cr \
    + (10.21 + 0.00835*T)*Mo \
    + (12.72 - 0.00325*T)*Ni

rhoF(T) = rhoF_(T, C, Cr, Mo, Ni)
rhoA(T) = rhoA_(T, C, Cr, Mo, Ni)
rhoL(T) = rhoL2(T, C, Cr, Mo, Ni)

print "rhoF = ", rhoF(TL)
print "rhoA = ", rhoA(TL)
print "rhoL = ", rhoL(TL)

range0(T, y) = (T > TS) ? 1/0 : y
range1(T, y) = (T > TS2) || (T < TS) ? 1/0 : y
range2(T, y) = (T > TA) || (T < TS2) ? 1/0 : y
range3(T, y) = (T > TL) || (T < TA) ? 1/0 : y
fcc(T, y) = (T > TA) ? 1/0 : y

array solute[4] = ['C', 'Cr', 'Mo', 'Ni']
array keys[4] = [ 'r t', 'r b', 'r b', 'r t']

array equilL[4] = [0.0003, 0.17, 0.025, 0.12]
array slopeL[4] = [-63800., 6272., 7877., -862.2]
array equilS[4] = [3.27e-5, 0.1747, 0.02881, 0.09136]
array slopeS[4] = [-674700., 5765., 6990., -1226.]

array equilL2[4] = [-5.4e-4, 0.1567, 0.01881, 0.1311]
array slopeL2[4] = [-12100., -1106., -2730., -3642.]
array equilS2[4] = [-4.74e-5, 0.1616, 0.02243, 0.09910]
array slopeS2[4] = [-127000., -1148., -2739., -5032.]

array equilL1[4] = [-5.4e-4, 0.1610, 0.02023, 0.13405] # Ni: 0.1402
array slopeL1[4] = [-12100., -1461., -3401., -8000.] # Ni: 5378
array equilS1[4] = [-4.74e-5, 0.1685, 0.02506, 0.10055] # Ni: 0.1050
array slopeS1[4] = [-127000., -1918., -4320., -8000.] # Ni: 10020

array equilS0[4] = [1.96e-4, 0.1790, 0.0297, 0.10269] # Ni: 0.1067
array slopeS0[4] = [1494000., -6527., -14790., -20000.] # Ni: 6214

do for [i=1:4] {
    j = i+2
    titleL = sprintf('Liquidus slope %.3g K/wt%', slopeL[i]/100)
    titleS = sprintf('Solidus slope %.3g K/wt%', slopeS[i]/100)
    titleP = "Piecewise approximation"
    eval('set key '.keys[i])
    plot 'liquid.txt' u ($1+DT):j w l, \
        'ferrite.txt' u ($1+DT):j w l, \
        equilL[i] + (x-TL)/slopeL[i] w l dt 2 lc 1 title titleL, \
        equilS[i] + (x-TL)/slopeS[i] w l dt 2 lc 2 title titleS, \
        range2(x, equilL2[i] + (x-TL)/slopeL2[i]) w l dt 3 lc 1 lw 2 notitle, \
        range2(x, equilS2[i] + (x-TL)/slopeS2[i]) w l dt 3 lc 2 lw 2 notitle, \
        range1(x, equilL1[i] + (x-TL)/slopeL1[i]) w l dt 3 lc 1 lw 2 notitle, \
        range1(x, equilS1[i] + (x-TL)/slopeS1[i]) w l dt 3 lc 2 lw 2 notitle, \
        range0(x, equilS0[i] + (x-TL)/slopeS0[i]) w l dt 3 lc 2 lw 2 title titleP, \
        'austenite.txt' u ($1+DT):(fcc($1+DT,column(j))) w l lc 5 #, \
#        '< paste ferrite.txt austenite.txt phases.txt' \
#            u ($1+DT):(column(j)*$15/($15+$16) + column(j+6)*$16/($15+$16)) w l dt 4 lc 4 \
#            title solute[i]."(bcc+fcc)"
}

unset multiplot

if (EPS) {
    set terminal postscript eps color font 'Helvetica,14' size 3.5in, 2in
    set output "2.eps"
} else {
    set term qt 1 position 900, 200
}

set ylabel "Phase fraction"
set key top center
plot for [i=2:4] 'phases.txt' u ($1+DT):i w l

