# gnuplot -c rep-seconds-flags.plt
# set terminal pdf enhanced color size 3.1482in,2.1318in font ",10"
set terminal pdf enhanced color size 3.1482in,2.1318in font ",10"
bench_type = "all"
y_index = 3
y_label = "seconds"
pdf_name = "rep-".bench_type."-seconds-flags.pdf"

set output pdf_name

load "common.plt"

array flags = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
array flags = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

impl = "vec4"

set margins 3, 1.75, 1.5, 3

datafile(i) = "flags-".i."-".bench_type.".dat"

firstrow = system('head -1 '.datafile(impl))

y_axis_pos = -0.01
set label "{/=12:Bold Sphere Trace (single precision) on Skylake 2.60GHz}" at character 0.01, screen 0.95

set xtics 100
set ytics 2

set label "Runtime [".y_label."] vs. number of shapes (mixed scene)" at character 0.01, screen 0.88

plot for [j=1:|flags|] datafile(impl) using 1:(column(3 * (flags[j] - 1) + y_index)) with linespoints linestyle j notitle

# vim:set ft=gnuplot:
