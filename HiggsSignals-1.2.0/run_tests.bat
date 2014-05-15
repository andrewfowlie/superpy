#!/bin/bash

startdir="$PWD"

# Checking for a eps viewer
if hash gv 2>/dev/null; then
	openeps=gv
elif hash open 2>/dev/null; then
	openeps=open
elif hash ggv 2>/dev/null; then
	openeps=ggv
elif hash eog 2>/dev/null; then
	openeps=eog
else
	echo 'Note: No eps viewer found. The .eps files will not be displayed.'
	openeps='#'
fi	

# Checking for a gnuplot
if hash gnuplot 2>/dev/null; then
	gnplt=gnuplot
else
	echo 'Note: Gnuplot not found. Will not run plotting scripts.'
	gnplt='#'
fi	



echo '*********************************************************************'
echo '*                   Running HiggsSignals test script                *'
echo '*********************************************************************'
echo 'cleaning HiggsSignals distribution...'
rm temp_error*.txt temp_output*.txt > /dev/null
make hyperclean 2>> temp_error.txt 1>> temp_output.txt
echo 'running configure script...'
./configure 2>> temp_error.txt 1>> temp_output.txt
echo 'make HiggsSignals...'
make HiggsSignals 2>> temp_error.txt 1>> temp_output.txt
if [ $(grep -c 'Error' temp_error.txt) == 0 ]; then
   echo 'Compilation successfully.'
else
   echo 'Compilation failed with'
   grep 'Error' temp_error.txt
   echo '---------------------------------------------------------------------'
   echo 'Here are some things you might want to check:'
   echo ' - Have you set the HiggsBounds path correctly in the configure script?'
   echo ' - Does the HiggsBounds library exist?'
   echo ' - Is it compiled with the same Fortran compiler?'   
   echo '---------------------------------------------------------------------'
   exit 0
fi   
echo '---------------------------------------------------------------------'
echo ' TEST 1: Run HiggsSignals command-line version on random test points:'
echo '---------------------------------------------------------------------'
echo './HiggsSignals latestresults peak 2 effC 3 1 example_data/random/HB_randomtest50points_'
./HiggsSignals latestresults peak 2 effC 3 1 example_data/random/HB_randomtest50points_ 2>> temp_error1.txt 1>> temp_output1.txt
# echo './HiggsSignals latestresults peak 2 part 3 1 example_data/random/HB_randomtest50points_'
# ./HiggsSignals latestresults peak 2 part 3 1 example_data/random/HB_randomtest50points_ 2>> temp_error1.txt 1>> temp_output1.txt
# echo './HiggsSignals latestresults peak 2 hadr 3 1 example_data/random/HB_randomtest50points_'
# ./HiggsSignals latestresults peak 2 hadr 3 1 example_data/random/HB_randomtest50points_ 2>> temp_error1.txt 1>> temp_output1.txt
# echo './HiggsSignals latestresults mass 2 effC 3 1 example_data/random/HB_randomtest50points_'
# ./HiggsSignals latestresults mass 2 effC 3 1 example_data/random/HB_randomtest50points_ 2>> temp_error1.txt 1>> temp_output1.txt
echo './HiggsSignals latestresults mass 2 part 3 1 example_data/random/HB_randomtest50points_'
./HiggsSignals latestresults mass 2 part 3 1 example_data/random/HB_randomtest50points_ 2>> temp_error1.txt 1>> temp_output1.txt
# echo './HiggsSignals latestresults mass 2 hadr 3 1 example_data/random/HB_randomtest50points_'
# ./HiggsSignals latestresults mass 2 hadr 3 1 example_data/random/HB_randomtest50points_ 2>> temp_error1.txt 1>> temp_output1.txt
# echo './HiggsSignals latestresults both 2 effC 3 1 example_data/random/HB_randomtest50points_'
# ./HiggsSignals latestresults both 2 effC 3 1 example_data/random/HB_randomtest50points_ 2>> temp_error1.txt 1>> temp_output1.txt
# echo './HiggsSignals latestresults both 2 part 3 1 example_data/random/HB_randomtest50points_'
/./HiggsSignals latestresults both 2 part 3 1 example_data/random/HB_randomtest50points_ 2>> temp_error1.txt 1>> temp_output1.txt
echo './HiggsSignals latestresults both 2 hadr 3 1 example_data/random/HB_randomtest50points_'
./HiggsSignals latestresults both 2 hadr 3 1 example_data/random/HB_randomtest50points_ 2>> temp_error1.txt 1>> temp_output1.txt

echo '---------------------------------------------------------------------'
echo 'There were' $(grep -c 'WARNING' temp_output.txt) 'warnings.'
echo 'There were' $(grep -c 'Interrupt' temp_error.txt) 'interrupts.'
echo 'There were' $(grep -c 'Error' temp_error.txt) 'errors.'
echo 'There were' $(grep -c 'stop' temp_error.txt) 'stops.'
echo '---------------------------------------------------------------------'
echo 'make HiggsSignals example programs...'
make HSexamples 2> temp_error.txt 1> temp_output.txt
if [ $(grep -c 'Error' temp_error.txt) == 0 ]; then
   echo 'Compilation successfully.'
else
   echo 'Compilation failed with'
   grep 'Error' temp_error.txt
fi
echo '---------------------------------------------------------------------'
echo ' TEST 2: Run HiggsSignals example programs:'
echo '---------------------------------------------------------------------'
echo -n '1) Running HSscaleUncertainties...'
cd example_programs
./HSscaleUncertainties 2 > ../temp_error2.txt 1 > ../temp_output2.txt
if [ $(grep -c 'Error' ../temp_error2.txt) == 0 ] || [ $(grep -c 'stop' ../temp_error2.txt) == 0 ]; then
   echo ' done.'
   echo -n '   Running gnuplot script...'
   cd results
   $gnplt plot_scaleUncertainties.gnu > /dev/null
   echo 'done (created example_programs/results/scaling_mu.eps).'
   $openeps scaling_dmu.eps 2>/dev/null &
   cd ..
else
   echo '   an error occured. Going to next example...'
fi
echo -n '2) Running HShadr...'
./HShadr 2 > ../temp_error2.txt 1 > ../temp_output2.txt
if [ $(grep -c 'Error' ../temp_error2.txt) == 0 ] || [ $(grep -c 'stop' ../temp_error2.txt) == 0 ]; then
   echo ' done.'
   echo -n '   Running gnuplot script...'
   cd results
   $gnplt plot_CSscaling.gnu > /dev/null
   echo 'done (created example_programs/results/CSscaling.eps).'
   $openeps CSscaling.eps 2>/dev/null &
   cd ..
else
   echo '   an error occured. Going to next example...'
fi
echo -n '3) Running HSeffC...'
./HSeffC 2 > ../temp_error2.txt 1 > ../temp_output2.txt
if [ $(grep -c 'Error' ../temp_error2.txt) == 0 ] || [ $(grep -c 'stop' ../temp_error2.txt) == 0 ]; then
   echo ' done.'
   echo -n '   Running gnuplot script...'
   cd results
   $gnplt plot_HSeffC.gnu > /dev/null
   echo 'done (created example_programs/results/Hgg_Hbb.eps).'
   $openeps Hgg_Hbb.eps 2>/dev/null &
   cd ..
else
   echo '   an error occured. Going to next example...'
fi
echo -n '4) Running HSwithSLHA on provided example SLHA file...'
./HSwithSLHA 1 ../example_data/SLHA/SLHA_FHexample.fh 2> ../temp_error2.txt 1> ../temp_output2.txt
if [ $(grep -c 'Error' ../temp_error2.txt) == 0 ] || [ $(grep -c 'stop' ../temp_error2.txt) == 0 ]; then
   echo ' done.'
else
   echo '   an error occured. Going to next example...'
fi
echo -n '5) Running HBandHSwithSLHA on provided example SLHA file...'
./HBandHSwithSLHA 1 ../example_data/SLHA/SLHA_FHexample.fh 2> ../temp_error2.txt 1> ../temp_output2.txt
if [ $(grep -c 'Error' ../temp_error2.txt) == 0 ] || [ $(grep -c 'stop' ../temp_error2.txt) == 0 ]; then
   echo ' done.'
else
   echo '   an error occured. Going to next example...'
fi
echo -n '6) Running HSwithToys...'
./HSwithToys 2 > ../temp_error2.txt 1 > ../temp_output2.txt
if [ $(grep -c 'Error' ../temp_error2.txt) == 0 ]; then
   echo ' done.'
else
   echo '   an error occured. Going to next example...'
fi
echo -n '7) Running HS_mass...'
./HS_mass 2 > ../temp_error2.txt 1 > ../temp_output2.txt
if [ $(grep -c 'Error' ../temp_error2.txt) == 0 ] || [ $(grep -c 'stop' ../temp_error2.txt) == 0 ]; then
   echo ' done.'
   echo -n '   Running gnuplot script...'
   cd results
   $gnplt plot_mh.gnu > /dev/null
   echo 'done (created example_programs/results/HS_mass.eps).'
   $openeps HS_mass.eps 2>/dev/null &
   cd ..
else
   echo '   an error occured. Going to next example...'
fi
echo -n '8) Running HS_efficiencies...'
./HS_efficiencies 2 > ../temp_error2.txt 1 > ../temp_output2.txt
if [ $(grep -c 'Error' ../temp_error2.txt) == 0 ] || [ $(grep -c 'stop' ../temp_error2.txt) == 0 ]; then
   echo ' done.'
   echo -n '   Running gnuplot script...'
   cd results
   $gnplt plot_efficiencies.gnu > /dev/null
   echo 'done (created example_programs/results/HS_efficiencies.eps).'
   $openeps HS_efficiencies.eps 2>/dev/null &
   cd ..
else
   echo '   an error occured.'
fi

echo '---------------------------------------------------------------------'
echo ' FINISHED WITH ALL TESTS. ENJOY!'
echo '---------------------------------------------------------------------'
rm -f results/tmp/*