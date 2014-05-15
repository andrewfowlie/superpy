#!/bin/bash

#first, need to compile and make the HiggsBounds library:
cd ..
make clean
./configure
make libHB

#now compile the program example-SM_vs_4thGen, which uses the HiggsBounds library 
#n.b. use the same compiler as you used to make the HiggsBounds library 
cd example_programs
gfortran example-SM_vs_4thGen.F -o example-SM_vs_4thGen -L.. -lHB

#run the program
./example-SM_vs_4thGen

echo ' **************************************************'
echo ' The output files are'
echo '   example-SM-results.dat'
echo '   example-4thGen-results.dat'
echo " and the key to the process numbers is in the file"
echo '   Key.dat'
echo ' **************************************************'
