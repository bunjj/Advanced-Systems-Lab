# gnuplot -c rep-flops.plt
set terminal pdf enhanced color size 3.1482in,2.1318in font ",10"
bench_type = "all"
y_index = 2
y_label = "GFlops"
pdf_name = "rep-".bench_type."-flops.pdf"

set output pdf_name

load "common.plt"
array impls = ["ref", "opt1", "opt3", "opt5", "vec4"]

set margins 5, 1.75, 1.5, 3

datafile = "-".bench_type.".dat"

firstrow = system('head -1 '.impls[1].datafile)

y_axis_pos = -0.01
set label "{/=12:Bold Sphere Trace (single precision) on Skylake 2.60GHz}" at character 0.01, screen 0.95

set xtics 100

set yrange [0:1000]
set ytics 200

set label "Operations [".y_label."] vs. number of shapes (mixed scene)" at character 0.01, screen 0.88

set style arrow 1 heads filled size 6,15

set label '{/:Italic Baseline}' at 180, 650 textcolor ls 1
set label '{/:Italic Inlining}' at 320, 640 textcolor ls 2
set label '{/:Italic Early Term.}' at 265, 450 textcolor ls 3
set label '{/:Italic MVMs}' at 320, 330 textcolor ls 4
set label '{/:Italic Vectorization}' at 220, 150 textcolor ls 5

set arrow from first 370, first 499 to first 370, first 824 as 1 ls 14
set label '1.8x' at 373, 643 textcolor ls 14
set arrow from first 370, first 339 to first 370, first 449 as 1 ls 14
set label '1.5x' at 373, 384 textcolor ls 14

plot for [j=1:|impls|] impls[j].datafile using 1:y_index with linespoints linestyle j notitle

# vim:set ft=gnuplot:
