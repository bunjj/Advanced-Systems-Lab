# gnuplot -c rep-flops-hist.plt
bench_type = "all"
pdf_name = "rep-".bench_type."-flops-hist.pdf"
title = "Flop Distributions in Different Implementations"
subtitle = "Operations [GFlops]"
load "rep-common.plt"

set output pdf_name

set lmargin 5
set rmargin 12

set key outside right top vertical Left reverse noenhanced autotitle columnhead nobox
set key invert samplen 4 spacing 1.23204 width 0 height 0 

datafile = "flops-".bench_type.".dat"

set style arrow 1 heads filled size 6,15

set style data histograms
set style histogram rowstacked

set yrange [0:1000]
set ytics 200
set style fill solid border -1
set boxwidth 0.75

set linetype 1 lc rgb "#282a2e"
set linetype 2 lc rgb "#2d2db9"
set linetype 3 lc rgb "#a03232"
set linetype 4 lc rgb "#008080"
set linetype 5 lc rgb "#65478f"
set linetype 6 lc rgb "#e69e00"
set linetype 7 lc rgb "#57b5e8"
set linetype 8 lc rgb "#e61f0f"
set linetype 9 lc rgb "cyan"
set linetype 10 lc rgb "sea-green"

plot for [j=3:12] datafile using (column(j) / 10**9):xtic(1) title column(j) fs pattern 3 lw 1

# vim:set ft=gnuplot:
