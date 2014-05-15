#!/bin/bash

#first, need to compile and make:
cd ..
make clean
./configure
make HBwithCPsuperH

#now run HBwithCPsuperH, using input contained in the file HBwithCPsuperH.input
cd example_programs
./HBwithCPsuperH < HBwithCPsuperH.input

echo ' *****************************************************'
echo ' The key to the process numbers is in the file Key.dat'
echo ' *****************************************************'
