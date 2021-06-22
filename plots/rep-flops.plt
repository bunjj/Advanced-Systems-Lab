# gnuplot -c rep-flops.plt
bench_type = "all"
y_index = 2
y_label = "GFlops"
pdf_name = "rep-".bench_type."-flops.pdf"
subtitle = "Operations [".y_label."] vs. number of shapes (mixed scene)"
load "rep-common.plt"

set output pdf_name

set lmargin 5

datafile = "-".bench_type.".dat"

set xtics 100

set yrange [0:1000]
set ytics 200

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
