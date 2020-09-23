#!/usr/bin/env gnuplot --persist -c
EPS=0

if (EPS) {
    set terminal postscript eps color font 'Helvetica,14'
    set output "tip.eps"
} else {
    set terminal qt size 1000, 600 font 'Helvetica,14'
}

set multiplot layout 2, 1
set fit errorvariables results

set xlabel "Time (ms)"

file = "log.solidificationFoam"
pulling_speed = 0.01
undercooling = 20
window = 1.5
array fit_range = [ 0.5, 0.95 ]

if (ARGC > 0) file = ARG1
if (ARGC > 1) pulling_speed = ARG2 + 0
if (ARGC > 2) undercooling = ARG3 + 0
if (ARGC > 3) window = ARG4 + 0

print "File: ", file
print "Pulling speed (m/s) = ", pulling_speed
print "Undercooling (K) = ", undercooling
print "Window = ", window

# Generate 4 columns: time, position, speed, undercooling
data = sprintf("<(awk '/^Time.*/{a=1e3*$3} /^Tip p/{b=$4} /^Tip s/{c=$4} /^Tip u/{print a,b,c,$4}' %s)", file)
xmax = system(sprintf("grep 'Tip s' %s | wc -l", file)) + 0
print "xmax = ", xmax

# Find maximum time when tips grow
set xrange [fit_range[1]*xmax:fit_range[2]*xmax]
pos(t) = C*t + D
fit pos(x) data u 0:1 via C,D
tmax = C*xmax + D
#plot data u 0:1 w l, C*x + D

win(t) = 1 + (window - 1)*exp(-t/tmax)
f(t, y, y0) = y < y0/win(t) ? 1/0 : (y > y0*win(t) ? 1/0 : y)

# Find steady-state speed and undercooling
fit A data u 0:3 via A
fit B data u 0:4 via B

do for [i=1:2] {
    set arrow from C*fit_range[i]*xmax + D, graph(0,0) \
                to C*fit_range[i]*xmax + D, graph(1,.75) nohead dt 5
}

set xrange [0:tmax]

set ylabel "Speed (m/s)"
plot data u 1:(f($1, $3, pulling_speed)) w l title 'Tip speed', \
    pulling_speed title "Pulling speed" lw 2, A t sprintf("Averaged value = %g +/- %.1g", A, A_err)

set ylabel "Undercooling (K)"
plot data u 1:(f($1, $4, undercooling)) w l title 'Tip undercooling', \
    B t sprintf("Averaged value = %g +/- %.1g", B, B_err)
