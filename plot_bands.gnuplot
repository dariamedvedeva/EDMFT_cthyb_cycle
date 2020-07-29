set encoding iso_8859
set terminal postscript eps enhanced color font 'Helvetica,20'
set output 'bands.eps'
unset xtics
unset ztics
unset key
set pm3d
set view map
unset border
unset colorbox
set xrange [ ] noextend
set yrange [ ] noextend
set cbrange [0:6]

#set palette rgb 21,22,23
set palette negative defined ( \
    0 '#D53E4F',\
    1 '#F46D43',\
    2 '#FDAE61',\
    3 '#FEE08B',\
    4 '#E6F598',\
    5 '#ABDDA4',\
    6 '#66C2A5',\
    7 '#3288BD' )

splot 'spectral.dat' u 1:2:(-$3/pi) w pm3d
