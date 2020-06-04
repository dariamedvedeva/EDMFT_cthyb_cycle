#!/bin/bash
gnuplot << 'EOF'

set terminal pdf
set output "plot.pdf"
#set multiplot layout 2,2

set xrange [0:40]
set yrange[:]

set title "Delta"
plot "Delta.dat" using 1:2 with lines lc 0 lw 2 dt '-' title 'Delta_{old} Re',\
    "Delta.dat" using 1:3 with lines lc 0 lw 2 title 'Delta_{old} Im',\
    "Delta_new.dat" using 1:2 with lines lc 7 lw 2 dt '-' title 'Delta_{new} Re',\
    "Delta_new.dat" using 1:3 with lines lc 7 lw 2 title 'Delta_{new} Im'
    
#set title "Lambda"
#plot "Phi_0.dat" using 1:2 with lines title 'Phi_0 Re',\
#    "Phi_0.dat" using 1:3 with lines title 'Phi_0 Im',\
#    "Lambda_new.dat" using 1:2 with lines title 'Lambda_{new} Re',\
#    "Lambda_new.dat" using 1:3 with lines title 'Lambda_{new} Im',\
#    "Lambda_new_smooth.dat" using 1:2 with lines title 'Lambda_{new}(smooth) Re',\
#    "Lambda_new_smooth.dat" using 1:3 with lines title 'Lambda_{new}(smooth) Im'
    
set title "1P GF"
plot "Gw.dat" using 1:2 with lines lc 0 lw 2 dt '-' title 'Gw Re',\
    "Gw.dat" using 1:3 with lines  lc 0 lw 2 title 'Gw Im',\
    "G_loc.dat" using 1:2 with lines  lc 7 lw 2 dt '-' title 'G_{loc} Re',\
    "G_loc.dat" using 1:3 with lines lc 7 lw 2 title 'G_{loc} Im'

#unset multiplot
set xrange [0:150]
set title "Self-energy"
plot "Sw.dat" using 1:2 with lines lc 0 title 'Sw Re',\
    "Sw.dat" using 1:3 with lines lc 0 title 'Sw Im',\
    "Sigma_smooth.dat" using 1:2 with lines lc 5 title 'Sw Re',\
    "Sigma_smooth.dat" using 1:3 with lines lc 5 title 'Sw Im',\
    
#set title "2P GF"
#plot "Xw.dat" using 1:2 with lines title 'Xw Re',\
#    "Xw.dat" using 1:3 with lines title 'Xw Im',\
#    "X_loc.dat" using 1:2 with lines title 'X_{loc} Re',\
#    "X_loc.dat" using 1:3 with lines title 'X_{loc} Im'

#unset multiplot
pause -1
'EOF'

