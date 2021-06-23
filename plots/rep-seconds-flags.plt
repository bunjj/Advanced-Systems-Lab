# gnuplot -c rep-seconds-flags.plt
bench_type = "all"
y_index = 3
pdf_name = "rep-".bench_type."-seconds-flags.pdf"
title = "Runtime with Different Compilers and Flags"
subtitle = "Runtime [seconds] vs. #shapes"

set output pdf_name

load "rep-common.plt"

array flags = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

impl = "vec4"

set lmargin 3

datafile(i) = "flags-".i."-".bench_type.".dat"

firstrow = system('head -1 '.datafile(impl))

set xtics 100
set ytics 2

plot for [j=1:|flags|] datafile(impl) using 1:(column(3 * (flags[j] - 1) + y_index)) with linespoints linestyle j notitle

# vim:set ft=gnuplot:
