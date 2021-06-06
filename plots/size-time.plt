# gnuplot -c size-time.plt <bench_type> <y-index> <y-label> <output-basename>

set terminal pdf enhanced color size 4.14in,3.38in font ",12"
bench_type = ARG1
y_index = ARG2 + 0
y_label = ARG3
base_name = ARG4
ps_name = base_name.".pdf"


set output ps_name

load "common.plt"
datafile = "-".bench_type.".dat"

set tmargin 2

set key left

firstrow = system('head -1 '.impls[1].datafile)

set offsets 0, 0, graph 0.3, graph 0

set xlabel word(firstrow, 1)
set label "[".y_label."]" at graph -0.071, graph 1.08

plot for [i=1:|impls|] impls[i].datafile using 1:y_index with linespoints linestyle i title impls[i]

# vim:set ft=gnuplot:
