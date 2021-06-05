set xtics axis nomirror out
set ytics scale 0

set tmargin 5

set border 1

set style line 1 lc rgb "#282a2e" lt 7 ps 0.6
set style line 2 lc rgb "#2d2db9" lt 7 ps 0.6
set style line 3 lc rgb "#a03232" lt 7 ps 0.6
set style line 4 lc rgb "#008080" lt 7 ps 0.6
set style line 5 lc rgb "#65478f" lt 7 ps 0.6
set style line 6 lc rgb "#e69e00" lt 7 ps 0.6
set style line 7 lc rgb "#57b5e8" lt 7 ps 0.6
set style line 8 lc rgb "#e61f0f" lt 7 ps 0.6

set style line 12 lc rgb "#ffffff" lt 1 lw 5
# For vertical lines
set style line 13 lc rgb "#000000" lw 3
set style line 14 lc rgb "#333333" lt 7 ps 0.6

set style arrow 1 heads filled size 8,15

set grid ytics ls 12

set obj rect from graph 0, graph 0 to graph 1, graph 1 fs noborder fc rgb "#e6e6e6" behind

array impls = ["ref", "opt0", "opt1", "opt3", "opt4", "opt5", "vec1", "vec2", "vec3", "vec4"]

# vim:set ft=gnuplot:
