#!/bin/bash
gnuplot << 'EOF'

set terminal pdf
set output "plot.pdf"


#  IN RANGE
set xrange [0:30]
#set yrange[-1.0:0.5]

set multiplot layout 1,2
set key font ",7"

set title "Delta"
plot "Delta.dat" using 1:2 with lines lc 0 lw 2 dt '-.-' title 'Delta Re',\
    "Delta.dat" using 1:3 with lines lc 0 lw 2 title 'Delta Im'

unset yrange
unset xrange
set title "Fourier transform"
plot "Delta_tau_ct_hyb.dat" using 1:2 with lines lc 0 lw 2 dt '-.-' title 'Delta_{tau} Re',\
    "Delta_tau_ct_hyb.dat" using 1:3 with lines lc 0 lw 2 title 'Delta_{tau} Im'
unset yrange

#
#set xrange [0:5]
#set title "Lambda"
#plot "Phi.dat" using 1:2 with lines lc 0 lw 2 dt '.'  title 'Phi_{old} Re',\
#    "Phi.dat" using 1:3 with lines lc 0 lw 2 title 'Phi_{old} Im',\
#    "Lambda_new.dat" using 1:2 with lines lc 2 lw 2 dt '-.-' title 'Lambda_{new} Re',\
#    "Lambda_new.dat" using 1:3 with lines lc 2 lw 2 title 'Lambda_{new} Im',\
#    "Lambda_new_smooth.dat" using 1:2 with lines lc 3 lw 2 dt '.' title 'Lambda_{new}(smooth) Re',\
#    "Lambda_new_smooth.dat" using 1:3 with lines lc 3 lw 2 title 'Lambda_{new}(smooth) Im'
#
#
#unset yrange
#unset xrange
#set title "Fourier for next iteration"
#plot "K_tau.dat" using 1:2 with lines lc 0 lw 2 dt '.'  title 'Phi_{old} Re',\
#    "K_tau.dat" using 1:3 with lines lc 0 lw 2 title 'Phi_{old} Im'

unset multiplot
pause -1
'EOF'

