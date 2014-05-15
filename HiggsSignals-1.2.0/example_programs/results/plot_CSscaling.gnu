set term postscript enhanced color eps

filename ='HShadr.dat'
filename_out = 'CSscaling.eps'

set xrange [0.:2.]
set yrange [0.:2.]
set dgrid3d 21,21,1

set pm3d map corners2color c1 clip4in

# set zrange [*:90.]
# set cbrange [*:90.]

# Find minimum
set output 'tmp/tmp.eps'
splot filename u ($2):($3):($6) notit w pm3d
min_z = GPVAL_DATA_Z_MIN
plot filename u ($2):($6 < min_z+0.0000001 ? $6 : 1/0) notit w p
min_pos_x = GPVAL_DATA_X_MIN
plot filename u ($3):($6 < min_z+0.0000001 ? $6 : 1/0) notit w p
min_pos_y = GPVAL_DATA_X_MIN

# print "      best-fit point:"
# print "      -------------- "
# print "      minimal Chi^2 Value = ", min_z
# print "      (mu_ggf, mu_VH) = (", min_pos_x,",", min_pos_y,")"

set contour
unset surface
set cntrparam bspline
unset clip
unset colorbox

set table 'tmp/1sigmacontour.dat'
set cntrparam levels discrete 2.2958
splot filename u ($2):($3):($6-min_z) notit w pm3d
unset table

set table 'tmp/2sigmacontour.dat'
set cntrparam levels discrete 5.99
splot filename u ($2):($3):($6-min_z) notit w pm3d
unset table

reset

set size 0.7,1

set pm3d map corners2color c1 clip4in
set palette rgbformulae 30,31,32

set xrange [0.:2.]
set yrange [0.:2.]
set dgrid3d 21,21,1

set xlabel '{/Symbol m}_{ggF+ttH} x B/B_{SM}'
set ylabel '{/Symbol m}_{VBF+VH} x B/B_{SM}' offset -1
set xtics -1.,1
set mxtics 5
set ytics -1.,1.
set mytics 5

set zrange [0:20]
set cbrange [0:20]

set output filename_out

set multiplot

set grid

set label 1 point ps 1.5 pt 3 lc rgb 'green' at min_pos_x, min_pos_y front
set label 2 point ps 1.5 pt 13 lc rgb '#F0E68C' at 1.0,1.0 front
set label 3 point ps 1.5 pt 12 lc rgb '#BB0000' at 1.0,1.0 front
set label 4 '{/Symbol D}{/Symbol c}^2' at 2.07,2.1 front

splot filename u ($2):($3):($6-min_z+0.0000001) notit w pm3d

set size 0.538,0.727
set origin 0.0837,0.15

unset xtics
unset ytics
unset xlabel
unset ylabel
unset clabel
unset label 1
unset label 2
unset label 3
unset label 4
unset surface
set cntrparam bspline
unset colorbox

#set key Left reverse spacing 1.5 box width 0 at 1.9, 1.9 opaque
set key Left reverse spacing 1.5 box width 0 at 1.9, 1.9

plot 'tmp/1sigmacontour.dat' u ($1):($2) w l lt 1 lw 4 lc rgb '#CCCCCC' title '68% C.L.',\
	 'tmp/2sigmacontour.dat' u ($1):($2) w l lt 2 lw 4 lc rgb '#AAAAAA' title '95% C.L.'
	 
unset multiplot	 

