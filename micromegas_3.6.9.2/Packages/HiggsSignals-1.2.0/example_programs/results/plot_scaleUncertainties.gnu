set term postscript enhanced color font ",20" eps

set output 'scaling_dmu.eps'

set xlabel 'Uncertainty scalefactor (in %)'
set ylabel '{/Symbol c}^2_{total}'
set grid

set xrange[50:100.]
set yrange[*:*]
plot 'scaling_dmu_exp.dat' u ($1)*100.0 : ($2) w l lt 1 lc 1 lw 4 title "only exp. uncertainty",\
     'scaling_dmu_th.dat' u ($1)*100.0 : ($2) w l lt 2 lc 2 lw 4 title "only th. uncertainty",\
     'scaling_dmu_both.dat' u ($1)*100.0 : ($2) w l lt 3 lc 3 lw 4 title "exp.+th. uncertainty"

