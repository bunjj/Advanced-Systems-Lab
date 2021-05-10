set xtics axis nomirror out
set ytics scale 0

set offsets 0, 0, graph 0.05, graph 0
set tmargin 5

set border 1

set style line 1 lc rgb "#008080" lt 7 ps 0.8
set style line 2 lc rgb "#2d2db9" lt 7 ps 0.8
set style line 3 lc rgb "#a03232" lt 7 ps 0.8
set style line 4 lc rgb "#282a2e" lt 7 ps 0.8
set style line 5 lc rgb "#65478f" lt 7 ps 0.8
set style line 6 lt 7 ps 0.8
set style line 12 lc rgb "#ffffff" lt 1 lw 5
# For vertical lines
set style line 13 lc rgb "#000000" lw 3
set grid ytics ls 12

set obj rect from graph 0, graph 0 to graph 1, graph 1 fs noborder fc rgb "#e6e6e6" behind

# vim:set ft=gnuplot: