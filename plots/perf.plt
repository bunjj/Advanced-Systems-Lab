# gnuplot -c perf.plt <bench_type>

set terminal pdf enhanced color size 4.14in,3.38in font ",12"
bench_type = ARG1
y_index = 2
pdf_name = "perf-".bench_type.".pdf"

set output pdf_name

load "common.plt"
datafile = "-".bench_type.".dat"

set tmargin 2

set key left

firstrow = system('head -1 '.impls[1].datafile)

y_axis_pos = 0

set label "[flops/cycle]" at graph y_axis_pos, graph 1.08

set yrange [0:5]

set arrow 102 from first 0, first 4 to first 400, 4 nohead
set label 102 '{/=10:Italic Peak Scalar (4 flops/cycle)}' at 270, 4.2

set xlabel word(firstrow, 1)

array impls = ["ref", "opt1"]
array lss = [1, 2]

plot for [i=1:|impls|] impls[i].datafile using 1:5 with linespoints linestyle lss[i] title impls[i]

unset label 100
unset label 101
unset label 102

set yrange [0:40]

set arrow from first 0, first 32 to first 400, 32 nohead
set label '{/=10:Italic Peak Vect. (32 flops/cycle)}' at 270, 33.5

array impls = ["opt5", "vec3"]
array lss = [4, 5]

set label 102 '{/=10:Italic Peak Scalar (4 flops/cycle)}' at 270, 6.5

plot for [i=1:|impls|] impls[i].datafile using 1:5 with linespoints linestyle lss[i] title impls[i]

# vim:set ft=gnuplot:
