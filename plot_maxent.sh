#!/bin/bash
gnuplot << 'EOF'

set terminal pdf
set output "plot_maxent.pdf"

#set multiplot layout 4,3

set title "in.out.avspec.dat"
plot "in.out.avspec.dat" u 1:2 w l, "in.out.avspec.dat" u 1:3 w l

set title "in.out.avspec_back.dat"
plot "in.out.avspec_back.dat" u 1:2 w l

set title "in.out.chi2.dat"
plot "in.out.chi2.dat" u 1:2 w l

set title "in.out.chispec.dat"
plot "in.out.chispec.dat" u 1:2 w l, "in.out.chispec.dat" u 1:3 w l

set title "in.out.chispec_back.dat"
plot "in.out.chispec_back.dat" u 1:2 w l

set title "in.out.fits.dat"
plot "in.out.fits.dat" u 1:2 w l, "in.out.fits.dat" u 1:3 w l, "in.out.fits.dat" u 1:4 w l

set title "in.out.maxspec.dat"
plot "in.out.maxspec.dat" u 1:2 w l

set title "in.out.maxspec_back.dat"
plot "in.out.maxspec_back.dat" u 1:2 w l

set title "in.out.prob.dat"
plot "in.out.prob.dat" u 1:2 w l

set title "in.out.spex.dat"
plot "in.out.spex.dat" u 1:2 w l

#unset multiplot

pause -1
'EOF'

