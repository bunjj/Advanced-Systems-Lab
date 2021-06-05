# gnuplot -c presentation.plt <bench_type> <y-index> <y-label> <output-basename>
set terminal pdf enhanced color size 4.14,3.38 font ",12"
bench_type = ARG1
y_index = ARG2 + 0
y_label = ARG3
base_name = ARG4
pdf_name = base_name.".pdf"

set output pdf_name

load "common.plt"
array impls = ["ref", "opt1", "opt3", "opt5", "vec2", "vec3"]

datafile = "-".bench_type.".dat"

set tmargin 2

firstrow = system('head -1 '.impls[1].datafile)

# set offsets 0, 0, graph 0.3, graph 0

set xtics 100

set yrange [0:400]
set ytics 100

# set logscale y 2
# set yrange [8:512]
# set ytics 2

# set xlabel word(firstrow, 1)
set label "[".y_label."]" at graph -0.075, graph 1.08

do for [i=1:|impls|] {
    plot for [j=1:i] impls[j].datafile using 1:y_index with linespoints linestyle j notitle
}

# vim:set ft=gnuplot:
