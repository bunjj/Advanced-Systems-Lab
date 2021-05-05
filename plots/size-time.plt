# This is the width of the text in my latex file
set terminal post enhanced color 
set output "size-time.ps"

datafile="size-time.dat"
firstrow = system('head -1 '.datafile)

set xtics axis nomirror out
set ytics scale 0

set offsets 0, 0, graph 0.02, graph 0
set tmargin 5

set xlabel "#Pixels"
set label '[seconds]' at graph -0.055, graph 1.08

set border 1

set style line 1 lc rgb "#008080" lt 7 ps 0.8
set style line 2 lc rgb "#2d2db9" lt 7 ps 0.8
set style line 3 lc rgb "#a03232" lt 7 ps 0.8
set style line 12 lc rgb "#ffffff" lt 1 lw 5
# For vertical lines
set style line 13 lc rgb "#000000" lw 3
set grid ytics ls 12

set obj rect from graph 0, graph 0 to graph 1, graph 1 fs noborder fc rgb "#e6e6e6" behind

plot datafile using 1:3 with linespoints ls 1 title "-O3"

system("ps2pdf size-time.ps")

# vim:set ft=gnuplot:
