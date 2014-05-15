#!/bin/bash

#first, need to compile and make:

cd ..	
make clean
./configure
make

cd example_programs

# now, pick command line options for HiggsBounds:

# specify the prefix for the input files and the number of neutral Higgs

prefix='../example_data/HB_randomtest50points_'
nH='3'
nHplus='1'

whichinput='part'   # whichinput can be 'part', 'hadr' or 'effC'
whichanalyses='LandH'   # whichanalyses can be 'LandH', 'onlyL', 'onlyH' or 'onlyP' 

# and run HiggsBounds:
 echo ../HiggsBounds $whichanalyses $whichinput $nH $nHplus $prefix
../HiggsBounds $whichanalyses $whichinput $nH $nHplus $prefix

echo ' **************************************************'
echo ' The output files are'
echo ' '"$prefix"HiggsBounds_results.dat
echo " and"
echo ' '"$prefix"Key.dat
echo ' **************************************************'
