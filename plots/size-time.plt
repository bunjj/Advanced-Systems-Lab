# This is the width of the text in my latex file
set terminal post enhanced color 
set output "size-time.ps"

datafile="bench.dat"
firstrow = system('head -1 '.datafile)

load "common.plt"

set offsets 0, 0, graph 0.3, graph 0

set label "{/=22:Bold Sphere Tracing on Haswell 3.30GHz}" at graph -0, graph 1.18

set xlabel "Width/Height [Pixels]"
set label '[seconds]' at graph -0.051, graph 1.08

array opts = [0, 1, 2, 3, 4]

plot for [i=1:|opts|] c = 3 * opts[i] + 3 datafile using 1:c with linespoints linestyle i title columnhead(c)

system("ps2pdf size-time.ps")

# vim:set ft=gnuplot:
