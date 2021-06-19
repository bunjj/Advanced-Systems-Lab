# gnuplot -c rep-perf.plt

set terminal pdf enhanced color size 3.1482in,2.1318in font ",10"
bench_type = "all"
pdf_name = "rep-".bench_type."-perf.pdf"

set output pdf_name

load "common.plt"
datafile = "-".bench_type.".dat"

firstrow = system('head -1 '.impls[1].datafile)

set label "{/=12:Bold Sphere Trace (single precision) on Skylake 2.60GHz}" at character 0.01, screen 0.95
set label "Performance [flops/cycle] vs number of shapes (mixed scene)" at character 0.01, graph 1.09

set margins 2, 1.75, 1.5, 3

set yrange [0:7]

set arrow 103 from first 0, first 4 to first 400, 4 nohead
set label 102 '{/=10:Italic Peak Scalar (4 flops/cycle)}' at 250, 4.3

array impls = ["ref", "opt1", "vec0"]
array lss = [1, 2, 3]

set label 100 '{/:Italic Baseline}' at 250, 1.4 textcolor ls lss[1]
set label 101 '{/:Italic Inlining}' at 250, 2.5 textcolor ls lss[2]
set label 105 '{/:Italic Inlining Vectorized}' at 250, 6.3 textcolor ls lss[3]

plot for [i=1:|impls|] impls[i].datafile using 1:5 with linespoints linestyle lss[i] notitle

unset label 105

set margins 3, 1.75, 1.5, 3.05

set yrange [0:12]
set ytics 2

array impls = ["opt2", "vec1"]
array lss = [4, 5]

set label 100 '{/:Italic Avoid sqrts}' at 250, 3 textcolor ls lss[1]
set label 101 '{/:Italic Vectorized}' at 250, 11 textcolor ls lss[2]
set label 102 '{/=10:Italic Peak Scalar (4 flops/cycle)}' at 250, 4.7

plot for [i=1:|impls|] impls[i].datafile using 1:5 with linespoints linestyle lss[i] notitle

array impls = ["opt3", "vec2"]
array lss = [4, 5]

set label 100 '{/:Italic Enclosing sphere}' at 250, 3.2 textcolor ls lss[1]
set label 101 '{/:Italic Vectorized}' at 250, 11 textcolor ls lss[2]

plot for [i=1:|impls|] impls[i].datafile using 1:5 with linespoints linestyle lss[i] notitle

array impls = ["opt5", "vec3"]
array lss = [4, 5]

set label 100 '{/:Italic Precompute MVMs}' at 250, 3 textcolor ls lss[1]
set label 101 '{/:Italic Vectorized}' at 250, 11 textcolor ls lss[2]

plot for [i=1:|impls|] impls[i].datafile using 1:5 with linespoints linestyle lss[i] notitle

# vim:set ft=gnuplot:
