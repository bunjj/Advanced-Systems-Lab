# gnuplot -c rep-seconds-vect.plt
bench_type = "all"
pdf_name = "rep-".bench_type."-seconds-vect.pdf"
subtitle = "Runtime [seconds] vs. number of shapes (mixed scene)"
load "rep-common.plt"

set output pdf_name

set lmargin 4

datafile = "-".bench_type.".dat"

set xtics 100

set style arrow 1 heads filled size 6,15

set arrow 100 from first 370, first 56 to first 370, first 155 as 1 ls 14
set label 101 '3.3x' at 373, 105.5 textcolor ls 14
set label 102 '{/:Italic Inlining}' at 250, 140 textcolor ls 14
set label 103 '{/:Italic Vectorized}' at 250, 60 textcolor ls 14

set ytics 40
plot "opt1".datafile using 1:3 with linespoints linestyle 1 notitle, \
"vec0".datafile u 1:3 with linespoints linestyle 3 notitle

set arrow 100 from first 370, first 36 to first 370, first 138 as 1 ls 14
set label 101 '5.0x' at 373, 87 textcolor ls 14
set label 102 '{/:Italic Avoid sqrts}' at 250, 130 textcolor ls 14
set label 103 '{/:Italic Vectorized}' at 250, 33 textcolor ls 14

set ytics 40
plot "opt2".datafile using 1:3 with linespoints linestyle 1 notitle, \
"vec1".datafile u 1:3 with linespoints linestyle 3 notitle


set arrow 100 from first 370, first 22 to first 370, first 79 as 1 ls 14
set label 101 '4.6x' at 373, 50.5 textcolor ls 14
set label 102 '{/:Italic Enclosing sphere}' at 200, 70 textcolor ls 14
set label 103 '{/:Italic Vectorized}' at 200, 27 textcolor ls 14

set ytics 20
plot "opt3".datafile using 1:3 with linespoints linestyle 1 notitle, \
"vec2".datafile u 1:3 with linespoints linestyle 3 notitle


set arrow 100 from first 370, first 14 to first 370, first 55 as 1 ls 14
set label 101 '4.7x' at 373, 34.5 textcolor ls 14
set label 102 '{/:Italic Precompute MVMs}' at 200, 45 textcolor ls 14
set label 103 '{/:Italic Vectorized}' at 200, 15 textcolor ls 14

set ytics 10
plot "opt5".datafile using 1:3 with linespoints linestyle 1 notitle, \
"vec3".datafile u 1:3 with linespoints linestyle 3 notitle

# vim:set ft=gnuplot:
