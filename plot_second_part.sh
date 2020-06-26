#!/bin/bash
gnuplot << 'EOF'

set terminal pdf
set output "plot.pdf"

set multiplot layout 1,2
#  IN RANGE
set xrange [0:5]
#set yrange[-1.0:0.5]

set title "1P GF"
plot "Gw.dat" using 1:2 with lines lc 0 lw 2 dt '-' title 'Gw Re',\
    "Gw.dat" using 1:3 with lines  lc 0 lw 2 title 'Gw Im',\
    "G_loc.dat" using 1:2 with lines  lc 7 lw 2 dt '-' title 'G_{loc} Re',\
    "G_loc.dat" using 1:3 with lines lc 7 lw 2 title 'G_{loc} Im'
    
set title "2P GF"
plot "Xw.dat" using 1:2 with lines lc 0 lw 2 dt '-.-' title 'Xw Re',\
    "Xw.dat" using 1:3 with lines lc 0 lw 2 title 'Xw Im',\
    "X_loc.dat" using 1:2 with lines lc 2 lw 2 dt '-' title 'X_{loc} Re',\
    "X_loc.dat" using 1:3 with lines lc 2 lw 2  title 'X_{loc} Im'
unset multiplot

set multiplot layout 2,2
set key font ",7"

set title "Delta"
plot "Delta.dat" using 1:2 with lines lc 0 lw 2 dt '-.-' title 'Delta_{old} Re',\
    "Delta.dat" using 1:3 with lines lc 0 lw 2 title 'Delta_{old} Im',\
    "Delta_new.dat" using 1:2 with lines lc 7 lw 2 dt '-.-' title 'Delta_{new} Re',\
    "Delta_new.dat" using 1:3 with lines lc 7 lw 2 title 'Delta_{new} Im',\
    "Delta_new_minimized.dat" using 1:2 with lines lc 6 lw 2 dt '-.-' title 'Delta_{min} Re',\
    "Delta_new_minimized.dat" using 1:3 with lines lc 6 lw 2 title 'Delta_{min} Im'

unset yrange
unset xrange
set title "Fourier for next iteration"
plot "Delta_new_minimized_tau_ct_hyb.dat" using 1:2 with lines lc 0 lw 2 dt '-.-' title 'Delta_{tau} Re',\
    "Delta_new_minimized_tau_ct_hyb.dat" using 1:3 with lines lc 0 lw 2 title 'Delta_{tau} Im'
unset yrange


set xrange [0:5]
set title "Lambda"
plot "Phi.dat" using 1:2 with lines lc 0 lw 2 dt '.'  title 'Phi_{old} Re',\
    "Phi.dat" using 1:3 with lines lc 0 lw 2 title 'Phi_{old} Im',\
    "Lambda_new.dat" using 1:2 with lines lc 2 lw 2 dt '-.-' title 'Lambda_{new} Re',\
    "Lambda_new.dat" using 1:3 with lines lc 2 lw 2 title 'Lambda_{new} Im',\
    "Lambda_new_smooth.dat" using 1:2 with lines lc 3 lw 2 dt '.' title 'Lambda_{new}(smooth) Re',\
    "Lambda_new_smooth.dat" using 1:3 with lines lc 3 lw 2 title 'Lambda_{new}(smooth) Im'


unset yrange
unset xrange
set title "Fourier for next iteration"
plot "K_tau.dat" using 1:2 with lines lc 0 lw 2 dt '.'  title 'Phi_{old} Re',\
    "K_tau.dat" using 1:3 with lines lc 0 lw 2 title 'Phi_{old} Im'

unset multiplot

set multiplot layout 2,2
#  WITHOUT RANGE
# # # # #
unset yrange
set xrange[0:45]
set title "Delta"
plot "Delta.dat" using 1:2 with lines lc 0 lw 2 dt '-.-' title 'Delta_{old} Re',\
    "Delta.dat" using 1:3 with lines lc 0 lw 2 title 'Delta_{old} Im',\
    "Delta_new.dat" using 1:2 with lines lc 7 lw 2 dt '-.-' title 'Delta_{new} Re',\
    "Delta_new.dat" using 1:3 with lines lc 7 lw 2 title 'Delta_{new} Im',\
    "Delta_new_minimized.dat" using 1:2 with lines lc 3 lw 2 dt '-.-' title 'Delta_{new} Re',\
    "Delta_new_minimized.dat" using 1:3 with lines lc 3 lw 2 title 'Delta_{new} Im'
unset yrange

set title "Lambda"
plot "Phi.dat" using 1:2 with lines lc 0 lw 2 dt '.'  title 'Phi_{old} Re',\
    "Phi.dat" using 1:3 with lines lc 0 lw 2 title 'Phi_{old} Im',\
    "Lambda_new.dat" using 1:2 with lines lc 2 lw 2 dt '-.-' title 'Lambda_{new} Re',\
    "Lambda_new.dat" using 1:3 with lines lc 2 lw 2 title 'Lambda_{new} Im',\
    "Lambda_new_smooth.dat" using 1:2 with lines lc 3 lw 2 dt '.' title 'Lambda_{new}(smooth) Re',\
    "Lambda_new_smooth.dat" using 1:3 with lines lc 3 lw 2 title 'Lambda_{new}(smooth) Im'
    
set title "1P GF"
plot "Gw.dat" using 1:2 with lines lc 0 lw 2 dt '-' title 'Gw Re',\
    "Gw.dat" using 1:3 with lines  lc 0 lw 2 title 'Gw Im',\
    "G_loc.dat" using 1:2 with lines  lc 7 lw 2 dt '-' title 'G_{loc} Re',\
    "G_loc.dat" using 1:3 with lines lc 7 lw 2 title 'G_{loc} Im'

   
set title "2P GF"
plot "Xw.dat" using 1:2 with lines lc 0 lw 2 dt '-.-' title 'Xw Re',\
    "Xw.dat" using 1:3 with lines lc 0 lw 2 title 'Xw Im',\
    "X_loc.dat" using 1:2 with lines lc 2 lw 2 dt '-' title 'X_{loc} Re',\
    "X_loc.dat" using 1:3 with lines lc 2 lw 2  title 'X_{loc} Im'
unset multiplot

set title "Self-energy"
plot "Sw.dat" using 1:2 with lines lc 0 title 'Sw Re',\
    "Sw.dat" using 1:3 with lines lc 0 title 'Sw Im'
#    "Sigma_smooth.dat" using 1:2 with lines lc 5 title 'Sw Re',\
#    "Sigma_smooth.dat" using 1:3 with lines lc 5 title 'Sw Im',\

pause -1
'EOF'

