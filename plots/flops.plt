# This is the width of the text in my latex file
set terminal post enhanced color 
set output "flops.ps"

datafile="bench.dat"

load "common.plt"

set label "{/=22:Bold Sphere Tracing Number of Flops}" at graph -0, graph 1.18

set xlabel "Width/Height [Pixels]"
set label '[Gflops]' at graph -0.051, graph 1.08

set label "{/:Bold Reference}" at 100, 15 textcolor ls 1

plot datafile using 1:($2 / 10 ** 9) with linespoints linestyle 1 notitle

system("ps2pdf flops.ps")

# vim:set ft=gnuplot:
