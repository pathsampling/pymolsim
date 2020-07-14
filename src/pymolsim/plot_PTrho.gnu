set lmargin 10
set bmargin 3
set rmargin 2
set tmargin 3

#set output 'W-BOP-Mirek-newmodelI-3x3x3.eps'
#set terminal postscript color

set multiplot

set title "Temperature"
set size 0.35,0.9
set origin 0.0,0.0
#set yrange [-0.5:3.0]
set ylabel 'temperature'
set xlabel 'time'
plot 'ensemble.dat' us 1:3 w lp lw 1 t 'T - average', \
     'ensemble.dat' us 1:2 w lp lw 1 t 'T - instantaneous'

set title "Pressure"
set origin 0.35,0.00
set ylabel 'pressure'
set xlabel 'time'
#set yrange [-0.6:0.6]
#set logscale y
#set format y '%3.0e'
plot 'ensemble.dat' us 1:5 w lp lw 1 t 'P - average', \
     'ensemble.dat' us 1:4 w lp lw 1 t 'P - instantaneous'



set title "Density"
set origin 0.7,0.0
set ylabel 'density'
set xlabel 'time'
#set yrange [-0.6:0.6]
#set logscale y
#set format y '%3.0e'
plot 'ensemble.dat' us 1:9 w lp lw 1 t 'rho - average', \
     'ensemble.dat' us 1:8 w lp lw 1 t 'rho - instantaneous'


#unset multiplot

#unset output
#set terminal x11
#pause -1 "Hit return to continue\n"

