# gnuplot -c pres-seconds.plt <bench_type>
set terminal pdf enhanced color size 4.14in,3.38in font ",12"
bench_type = ARG1
y_index = 3
y_label = "seconds"
pdf_name = "pres-".bench_type."-seconds.pdf"

set output pdf_name

load "common.plt"
array impls = ["ref", "opt1", "opt3", "opt5", "vec3"]

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

set xlabel word(firstrow, 1)
set label "[".y_label."]" at graph -0.075, graph 1.08

array label_eval = [ \
    "set label '{/:Italic Baseline}' at 250, 320 textcolor ls 1", \
    "set label '{/:Italic Inlining}' at 270, 150 textcolor ls 2", \
    "set label '{/:Italic Early Term.}' at 270, 85 textcolor ls 3", \
    "set label '{/*.8:Italic MVMs}' at 295, 56 textcolor ls 4", \
    "set label '{/:Italic Vectorization}' at 260, 25 textcolor ls 5" \
]

set style arrow 1 heads filled size 8,15

array arrows = [ \
    "set arrow from first 370, first 165 to first 370, first 365 as 1 ls 14", \
    "set arrow from first 370, first 89 to first 370, first 151 as 1 ls 14", \
    "set arrow from first 370, first 63 to first 370, first 75 as 1 ls 14", \
    "set arrow from first 370, first 18 to first 370, first 49 as 1 ls 14" \
]

array arrow_labels = [ \
    "set label '2.4x' at 373, 265 textcolor ls 14", \
    "set label '1.9x' at 373, 120 textcolor ls 14", \
    "set label '1.5x' at 373, 69 textcolor ls 14", \
    "set label '4.7x' at 373, 33 textcolor ls 14" \
]

do for [i=1:|impls|] {
    eval label_eval[i];

    if (i > 1) {
        eval arrows[i - 1];
        eval arrow_labels[i - 1];
    }
    plot for [j=1:i] impls[j].datafile using 1:y_index with linespoints linestyle j notitle
}

# vim:set ft=gnuplot:
