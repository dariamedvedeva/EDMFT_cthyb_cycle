#!/bin/bash
gnuplot << 'EOF'

set terminal pdf
set output "plot.pdf"
set multiplot layout 2,2

set xrange [0:5]
set yrange[: 0.25]

set title "Delta"
plot "Delta_0.dat" using 1:2 with lines title 'Delta_0 Re',\
    "Delta_0.dat" using 1:3 with lines title 'Delta_0 Im',\
    "Delta_new.dat" using 1:2 with lines title 'Delta_{new} Re',\
    "Delta_new.dat" using 1:3 with lines title 'Delta_{new} Im'

set title "Lambda"
plot "Phi_0.dat" using 1:2 with lines title 'Phi_0 Re',\
    "Phi_0.dat" using 1:3 with lines title 'Phi_0 Im',\
    "Lambda_new.dat" using 1:2 with lines title 'Lambda_{new} Re',\
    "Lambda_new.dat" using 1:3 with lines title 'Lambda_{new} Im',\
    "Lambda_new_smooth.dat" using 1:2 with lines title 'Lambda_{new}(smooth) Re',\
    "Lambda_new_smooth.dat" using 1:3 with lines title 'Lambda_{new}(smooth) Im'
    
set title "1P GF"
plot "Gw.dat" using 1:2 with lines title 'Gw Re',\
    "Gw.dat" using 1:3 with lines title 'Gw Im',\
    "G_loc.dat" using 1:2 with lines title 'G_{loc} Re',\
    "G_loc.dat" using 1:3 with lines title 'G_{loc} Im',\
    "Sw.dat" using 1:2 with lines title 'Sw Re',\
    "Sw.dat" using 1:3 with lines title 'Sw Im'
    
set title "2P GF"
plot "Xw.dat" using 1:2 with lines title 'Xw Re',\
    "Xw.dat" using 1:3 with lines title 'Xw Im',\
    "X_loc.dat" using 1:2 with lines title 'X_{loc} Re',\
    "X_loc.dat" using 1:3 with lines title 'X_{loc} Im'

unset multiplot
pause -1
'EOF'

