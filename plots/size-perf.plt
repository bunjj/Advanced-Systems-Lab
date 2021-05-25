set terminal post enhanced color 
set output "size-perf.ps"

datafile="bench.dat"
firstrow = system('head -1 '.datafile)

load "common.plt"

set offsets 0, 0, graph 0.2, graph 0

set label "{/=22:Bold Sphere Tracing on Haswell 3.30GHz}" at graph -0, graph 1.18

set xlabel "Width/Height [Pixels]"
set label '[flops/cycle]' at graph -0.071, graph 1.08

array opts = [0, 1, 2, 3, 4]

plot for [i=1:|opts|] c = 3 * opts[i] + 5 datafile using 1:c with linespoints linestyle i title columnhead(c)

system("ps2pdf size-perf.ps")

# vim:set ft=gnuplot:
