# gnuplot -c rep-seconds.plt
bench_type = "all"
y_index = 3
pdf_name = "rep-".bench_type."-seconds.pdf"
subtitle = "Runtime [seconds] vs. number of shapes (mixed scene)"

set output pdf_name

load "rep-common.plt"

set lmargin 4

datafile = "-".bench_type.".dat"

set xtics 100

set yrange [0:405]
set ytics 100

set style arrow 1 heads filled size 6,15

set label '{/:Italic Baseline}' at 250, 320 textcolor ls 1
set label '{/:Italic Inlining}' at 270, 150 textcolor ls 2
set label '{/:Italic Early Term.}' at 270, 87 textcolor ls 3
set label '{/*.7:Italic MVMs}' at 295, 56 textcolor ls 4
set label '{/:Italic Vectorization}' at 260, 28 textcolor ls 5

set arrow from first 370, first 172 to first 370, first 394 as 1 ls 14
set arrow from first 370, first 92 to first 370, first 152 as 1 ls 14
set arrow from first 370, first 67 to first 370, first 72 as 1 ls 14
set arrow from first 370, first 22 to first 370, first 47 as 1 ls 14

set label '2.5x' at 373, 283 textcolor ls 14
set label '2.0x' at 373, 122 textcolor ls 14
set label '1.5x' at 373, 69.5 textcolor ls 14
set label '4.7x' at 373, 34.5 textcolor ls 14

plot for [j=1:|impls|] impls[j].datafile using 1:y_index with linespoints linestyle j notitle

# vim:set ft=gnuplot:
