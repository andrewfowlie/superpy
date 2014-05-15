#!/bin/bash

#This script is part of HiggsBounds
# -KW

#It is called by configure. It creates the module store_pathname,
#which contains the path to the HiggsBounds package

#startdir="$PWD"
startdir=`pwd -P`
nameofpath=$startdir

teststring=$nameofpath

teststringlength=${#teststring}
maxlength=100
pathnamelength=`expr $teststringlength + 1` #because we're adding another '/' to the end

newstring=$teststring

echo '!******************************************************************'
echo 'module store_pathname'
echo '!******************************************************************'
echo ' implicit none'
echo ''
echo ' integer,parameter:: pathname_length= '$pathnamelength
echo ' character(len=pathname_length),parameter :: pathname= &'
echo '     &     "'${newstring:0:maxlength}'" // &'

newstring=${newstring:maxlength} 

while [ "${#newstring}" -gt "$maxlength" ]
do
 echo '     &     "'${newstring:0:maxlength}'" // &'
 newstring=${newstring:maxlength} 
done 

echo '     &     "'$newstring'/"'

echo ''
echo 'end module store_pathname'
echo '!******************************************************************'
