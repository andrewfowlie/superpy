/*
 Copyright (C) 1997, Alexander Pukhov 
*/

#include<stdlib.h>
#include<stdio.h>
#include<usrfun.h>

/*                FOR USER 

The given file contains a dummy version of 
                usrFF
function which should be  replaced by other ones written by the user.
We expect that user version will be attached to calchep numerical code
as  via 'Libraries' model file, while the current dummy version will 
be kept unchanged here. 

The n_in/n_out parameters present numbers in incoming/outgoing particles 
correspondiongly.
 
pvect  present momenta of particles 

             q[k]=pvect[4*(I-1)+k]  k=0,1,2,3 - momenta of I^{th} particle; 
                  1<=I<=nin_int          - incoming particles;
             nin_int<I<=nin_int+nout_int - outgoing particles;
             Energy of all particle are positive  pvect[4*(I-1)]>0;
             Axis of collision k=3.    
pnames array contains names on paricles in reaction. For example, 
pname[0] - in the name of first particle with energy pvect[0]. 
 
Below we present tools which can be used for  usrfun realization:        
*/

                                                                
double usrFF(int n_in, int n_out, double * pvect,char**pnames,int*pdg) { return 1.; }
