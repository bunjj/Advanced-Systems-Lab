# gnuplot speedup.plt 
set terminal pdf enhanced color size 4.14in,3.38in font ",12"
set output "speedup.pdf"

load "common.plt"
array impls = ["ref", "opt1", "opt3", "opt5", "vec3"]
array bench_types = ["size", "all", "box", "cone", "octahedron", "sphere", "torus"]

set tmargin 2


# set offsets 0, 0, graph 0.3, graph 0
set key left

set label "[speedup]" at graph 0, graph 1.08

do for [i=1:|bench_types|] {
    datafile = "-".bench_types[i].".dat"

    firstrow = system('head -1 '.impls[1].datafile)
    set xlabel word(firstrow, 1)

    set table $reference
    plot impls[1].datafile using 3 with table
    unset table


    plot for [j=2:|impls|] impls[j].datafile using 1:($reference[$0]/column(3)) with linespoints linestyle j title impls[j]
}

# vim:set ft=gnuplot:
