#!/usr/bin/env -S gnuplot --persist -c

EPS=1

if (EPS) {
    set terminal postscript eps color font 'Helvetica,14'
    set output "tip.eps"
} else {
    set terminal qt size 1000, 600 font 'Helvetica,14'
}

set multiplot layout 2, 1
set fit errorvariables results

set xlabel "Time (ms)"
time_factor = 1e3

log_file = "log.solidificationFoam"
array fit_range = [ 0.7, 1 ]

window = 2
if (ARGC > 0) window = ARG1 + 0

# Read problem properties
problem_file = "constant/problemProperties"
undercooling = system(sprintf("awk '/^undercooling/{print $2}' %s | sed 's/;//'", problem_file)) + 0
G = system(sprintf("awk '/^tempGradient/{print $2}' %s | sed 's/;//'", problem_file)) + 0
dotT = system(sprintf("awk '/^coolingRate/{print $2}' %s | sed 's/;//'", problem_file)) + 0
Vp = dotT/G

print "Initial undercooling (K) = ", undercooling
print "Temperature gradient (K/m) = ", G
print "Cooling rate (K/s) = ", dotT
print "Pulling speed (m/s) = ", Vp
print "Window = ", window

# Generate 4 columns: time, position, speed, undercooling
data = sprintf("<(awk '/^Time.*/{a=%f*$3} /^Tip p/{b=$4} /^Tip s/{c=$4} /^Tip u/{print a,b,c,$4}' %s)", time_factor, log_file)
xmax = system(sprintf("grep 'Tip s' %s | wc -l", log_file)) + 0
print "Number of points = ", xmax

# Find maximum time when tips grow
pos(t) = t < tmax ? C*t + D : C*tmax + D
stats data using 1 name "time" nooutput
tmax = time_max/2
fit pos(x) data u 1:2 via C, D, tmax
#plot data u 1:2 w lp, pos(x)

win(t) = 1 + (window - 1)*exp(-t/tmax)
f(t, y, y0) = y < y0/win(t) ? 1/0 : (y > y0*win(t) ? 1/0 : y)
set xrange [fit_range[1]*tmax:fit_range[2]*tmax]

# Find steady-state speed and undercooling
fit A data u 1:3 via A

f_undercooling(x) = B
fit B data u 1:4 via B

if (B_err/B > 1e-4) {
    set xrange [0.1*tmax:tmax]
    Q = 1
    f_undercooling(x) = B + P*exp(-x/Q**2)
    fit f_undercooling(x) data u 1:4 via B, P, Q
}

do for [i=1:2] {
    set arrow from fit_range[i]*tmax, graph(0,0) \
                to fit_range[i]*tmax, graph(1,.75) nohead dt 5
}

set xrange [0:tmax]

set ylabel "Speed (m/s)"
plot data u 1:(f($1, $3, Vp)) w l title 'Tip speed', \
    Vp title "Pulling speed" lw 2, A t sprintf("Averaged value = %g +/- %.1g", A, A_err)

set ylabel "Undercooling (K)"
plot data u 1:(f($1, $4, undercooling)) w l title 'Tip undercooling', \
    f_undercooling(x) t sprintf("Averaged value = %g +/- %.1g", B, B_err)
