set terminal post enhanced color 
bench_type = ARG1
y_index = ARG2 + 0
y_label = ARG3

base_name = "size-time-".bench_type."-".y_label
ps_name = base_name.".ps"

set output ps_name

load "common.plt"
datafile = "-".bench_type.".dat"

firstrow = system('head -1 '.impls[1].datafile)

set offsets 0, 0, graph 0.3, graph 0

set label "{/=22:Bold Sphere Tracing on Haswell 3.30GHz}" at graph -0, graph 1.18

set xlabel word(firstrow, 1)
set label "[".y_label."]" at graph -0.051, graph 1.08

plot for [i=1:|impls|] impls[i].datafile using 1:y_index with linespoints linestyle i title impls[i]

system("ps2pdf ".ps_name)

# vim:set ft=gnuplot:
