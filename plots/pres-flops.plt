# gnuplot -c pres-flops.plt <bench_type>
set terminal pdf enhanced color size 4.14in,3.38in font ",12"
bench_type = ARG1
y_index = 2
y_label = "GFlops"
pdf_name = "pres-".bench_type."-flops.pdf"

set output pdf_name

load "common.plt"
array impls = ["ref", "opt1", "opt3", "opt5", "vec3"]

datafile = "-".bench_type.".dat"

set tmargin 2

firstrow = system('head -1 '.impls[1].datafile)

# set offsets 0, 0, graph 0.3, graph 0

set xtics 100

# set yrange [0:400]
# set ytics 100

# set logscale y 2
# set yrange [8:512]
# set ytics 2

set xlabel word(firstrow, 1)
set label "[".y_label."]" at graph -0.075, graph 1.08

array label_eval = [ \
    "set label '{/:Italic Baseline}' at 180, 550 textcolor ls 1", \
    "set label '{/:Italic Inlining}' at 300, 640 textcolor ls 2", \
    "set label '{/:Italic Distance Functions}' at 220, 450 textcolor ls 3", \
    "set label '{/:Italic Precompute MVM}' at 300, 350 textcolor ls 4", \
    "set label '{/:Italic Vectorization}' at 220, 150 textcolor ls 5" \
]

do for [i=1:|impls|] {
    eval label_eval[i];
    plot for [j=1:i] impls[j].datafile using 1:y_index with linespoints linestyle j notitle
}

# vim:set ft=gnuplot:
