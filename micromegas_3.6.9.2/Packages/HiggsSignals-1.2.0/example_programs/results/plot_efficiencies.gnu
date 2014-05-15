set term postscript enhanced color font ",20" eps 

set output "HS_efficiencies.eps"

at(file, row, col) = system( sprintf("awk -v row=%d -v col=%d 'NR == row {print $col}' %s", row, col, file) )

file="HS_efficiencies_kv0p5.dat" ; row=41

# print "Reading in values with scale factor ",at(file,row,2)

sf0p5eff1 = at("HS_efficiencies_kv0p5.dat",row,3)
sf1p0eff1 = at("HS_efficiencies_kv1p0.dat",row,3)
sf1p5eff1 = at("HS_efficiencies_kv1p5.dat",row,3)
sf2p0eff1 = at("HS_efficiencies_kv2p0.dat",row,3)

set xrange [0.:2.5]

set xlabel "{/Symbol z} = {/Symbol e}^{model}/{/Symbol e}^{SM} (VBF, WH, ZH)"
set ylabel "{/Symbol c}^2 / {/Symbol c}^2_{({/Symbol z}=1)} -1"

set grid

set key Left reverse spacing 1.5 box width 0 at screen 0.4,0.8 opaque

set label 1 'Coupling scale factor {/Symbol k}_V varied.' at screen 0.2,0.9 font ',18'  #at 2.,-0.25  font ',18'
set label 2 'Others as in SM.' at screen 0.2,0.85 font ',18'

plot "HS_efficiencies_kv0p5.dat" u ($2):(($3)/sf0p5eff1 -1.0) w l lt 1 lc 1 lw 4 title "{/Symbol k}_V^2 = 0.5",\
      "HS_efficiencies_kv1p0.dat" u ($2):(($3)/sf1p0eff1 -1.0) w l lt 2 lc 2 lw 4 title "{/Symbol k}_V^2 = 1.0",\
      "HS_efficiencies_kv1p5.dat" u ($2):(($3)/sf1p5eff1 -1.0) w l lt 3 lc 3 lw 4 title "{/Symbol k}_V^2 = 1.5",\
      "HS_efficiencies_kv2p0.dat" u ($2):(($3)/sf2p0eff1 -1.0) w l lt 4 lc 4 lw 4 title "{/Symbol k}_V^2 = 2.0"