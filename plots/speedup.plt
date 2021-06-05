# gnuplot -c speedup.plt <bench_type> <y-index> <output-basename>
set terminal pngcairo enhanced color
bench_type = ARG1
y_index = ARG2 + 0
base_name = ARG3
pdf_name = base_name.".png"

set output pdf_name

load "common.plt"
array impls = ["ref", "opt1", "opt3", "opt5", "vec2", "vec3"]

datafile = "-".bench_type.".dat"

set tmargin 2

firstrow = system('head -1 '.impls[1].datafile)

set table $reference
plot impls[1].datafile using 3 with table
unset table

# set offsets 0, 0, graph 0.3, graph 0
set key outside

set xtics 100

set xlabel word(firstrow, 1)
set label "[speedup]" at graph -0.065, graph 1.08

plot for [i=2:|impls|] impls[i].datafile using 1:($reference[$0]/column(y_index)) with linespoints linestyle i title impls[i]

# vim:set ft=gnuplot:
