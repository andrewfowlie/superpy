/*
 Copyright (C) 1997, Alexander Pukhov 
*/

#include<stdlib.h>
#include<stdio.h>
#include<usrfun.h>


/*                FOR USER 

The given file contains a dummy version of 
                usrfun
function which should be  replaced by other ones written by the user.
We expect that user version will be attached to calchep numerical code
as  via 'Libraries' model file, while the current dummy version will 
be kept unchanged here. 
  If one uses U<txt> for  n_calchep  cuts and histograms then usrfun(<txt>) 
will be called. Say 'Uabs'  corresponds to usrfun("abs",pvect).

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
 
 
Below we present tools which can be used for usrfun.       
*/


                                                                
double usrfun(char * name,int n_in, int n_out, double * pvect,char**p_names,int*pdg)
{   
   fprintf(stdout," usrfun(char* name) is  called with parameter %s\n"
                  "But  is not defined!\n",name);
   exit(54);
   return 0.;
}
