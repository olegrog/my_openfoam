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
rhoL_(T, C, Cr, Mo, Ni) = 8319.49 - 0.835*T \
    + (-83.19 + 0.00835*T)*C \
    + (-14.77 + 0.00535*T)*Cr \
    + (10.21 + 0.00835*T)*Mo \
    + (12.72 - 0.00325*T)*Ni

rhoF(T) = rhoF_(T, C, Cr, Mo, Ni)
rhoA(T) = rhoA_(T, C, Cr, Mo, Ni)
rhoL(T) = rhoL_(T, C, Cr, Mo, Ni)

print "rhoF = ", rhoF(TL)
print "rhoA = ", rhoA(TL)
print "rhoL = ", rhoL(TL)

fcc(T, y) = T > TA ? 1/0 : y

array solute[4] = ['C', 'Cr', 'Mo', 'Ni']
array equilL[4] = [0.0003, 0.17, 0.025, 0.12]
array slopeL[4] = [-63803, 6271.86, 7877.23, -862.19]
array equilS[4] = [3.27e-5, 0.17471, 0.02881, 0.09136]
array slopeS[4] = [-674731, 5764.81, 6990.48, -1225.86]

array equilL_[4] = [-0.00054, 0.15819, 0.01932, 0.13433]
array slopeL_[4] = [-12106.98, -1279.01, -3067.18, -55294.12]
array equilS_[4] = [-4.74e-05, 0.16404, 0.02337, 0.10119]
array slopeS_[4] = [-127128, -1478.3, -3440, -30790]

do for [i=1:4] {
    j = i+2
    titleL = sprintf('Liquidus slope %.3g K/wt%', slopeL[i]/100)
    titleS = sprintf('Solidus slope %.3g K/wt%', slopeS[i]/100)
    plot 'liquid.txt' u ($1+DT):j w l, \
        'ferrite.txt' u ($1+DT):j w l, \
        equilL[i] + (x-TL)/slopeL[i] w l dt 2 lc 1 title titleL, \
        equilS[i] + (x-TL)/slopeS[i] w l dt 2 lc 2 title titleS, \
        fcc(x, equilL_[i] + (x-TL)/slopeL_[i]) w l dt 3 lc 1 notitle, \
        fcc(x, equilS_[i] + (x-TL)/slopeS_[i]) w l dt 3 lc 2 notitle, \
        'austenite.txt' u ($1+DT):(fcc($1+DT,column(j))) w l, \
        '< paste ferrite.txt austenite.txt phases.txt' \
            u ($1+DT):(column(j)*$15/($15+$16) + column(j+6)*$16/($15+$16)) w l dt 4 lc 4 \
            title solute[i]."(bcc+fcc)"
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

