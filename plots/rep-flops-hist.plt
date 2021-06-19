# gnuplot -c rep-flops-hist.plt
set terminal pdf enhanced color size 3.1482in,2.1318in font ",10"
bench_type = "all"
pdf_name = "rep-".bench_type."-flops-hist.pdf"

set output pdf_name

load "common.plt"
array impls = ["ref", "opt1", "opt3", "opt5", "vec4"]

set margins 5.5, 12.75, 1.5, 3.02

set key outside right top vertical Left reverse noenhanced autotitle columnhead nobox
set key invert samplen 4 spacing 1 width 0 height 0 

datafile = "flops-".bench_type.".dat"

set label "{/=12:Bold Sphere Trace (single precision) on Skylake 2.60GHz}" at character 0.01, screen 0.95

set label "Flop distribution" at character 0.01, screen 0.89

set style arrow 1 heads filled size 6,15

set style data histograms
set style histogram rowstacked

set yrange [0:100]
set format y "%g%%"
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

plot for [j=3:12] datafile using (100. * column(j) / $2):xtic(1) title column(j) fs pattern 3 lw 1

# vim:set ft=gnuplot:
