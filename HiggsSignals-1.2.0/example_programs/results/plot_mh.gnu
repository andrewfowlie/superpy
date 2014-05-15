set term postscript enhanced color font ",24"

set xrange [110.:140.]
set yrange [50:150]
set yrange [*:*]
set y2range [0:*]


set xtics 110.,5.
set mxtics 5
set mytics 5
set ytics nomirror
set y2tics 0,5 nomirror

set grid front

set xlabel 'm_{H} [GeV]'
set ylabel '{/Symbol c}^2'
set y2label 'Number of assignments'

# set key right font ',16' at 140,74 samplen 2 
set key bottom right font ',16' samplen 2 

file1 = "HS_mass_pdf1.dat"
file2 = "HS_mass_pdf2.dat"
file3 = "HS_mass_pdf3.dat"
fileout = "HS_mass.eps"

set output fileout
plot file1 u ($1):($7) axes x1y2 w l lt 1 lc rgb "#FFBBBB" lw 4 notit,\
     file2 u ($1):($7) axes x1y2 w l lt 2 lc rgb "#BBFFBB" lw 4 notit,\
     file3 u ($1):($7) axes x1y2 w l lt 3 lc rgb "#BBBBFF" lw 4 notit,\
     file1 u ($1):($6) w l lt 1 lc 1 lw 6 title "box",\
     file2 u ($1):($6) w l lt 2 lc 2 lw 6 title "Gaussian",\
     file3 u ($1):($6) w l lt 3 lc 3 lw 6 title "box+Gaussian"
