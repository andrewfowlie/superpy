#!/bin/bash

#first, need to compile and make:

cd ..
	
make clean
./configure
make HBwithFH

cd example_programs

./HBwithFH

echo ' *****************************************************'
echo ' The key to the process numbers is in the file Key.dat'
echo ' *****************************************************'
