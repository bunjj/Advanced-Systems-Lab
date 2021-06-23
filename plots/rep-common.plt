set terminal pdf enhanced color size 3.1482in,2.1318in font ",10"
load "common.plt"
array impls = ["ref", "opt1", "opt3", "opt5", "vec4"]

set rmargin 1.75
set bmargin 1.5
set tmargin 3

set label "{/=12:Bold ".title."}" at character 0.01, screen 0.95

set label subtitle at character 0.01, screen 0.88

# vim:set ft=gnuplot:
