#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<string.h>

extern double usrfun(char*name, int nIn, int nOut,  double * pvect,char**pName,int*pCode);
/*     
Any time we use  U<name>  to specify cuts and distribution, CalcHEP calls
function 'usrfun' to calculate corresponding value. CalcHEP code has 
a dummy version of this function stops calculation. The dummy version is 
replaced on the  user one if its code is passed to CalcHEP linker via 
'Libraries' model item.  One can use CALCHEP and  WORK  environment 
variables  to specify path to the code. These variables  are automatically 
defined in calchep and calchep_batch  scripts. Also one can use any other 
environment variables defined separately. 

 Parameters of usrFF:
    name - identifier of constructed physical variable. For instance, if 
       you call Uabs  in cuts/distibutions then  string "abs" is passed to 
       usrfun as the 'name' argument.  
    nIn - number of incoming particles;
    nOut- number of outgoing particles;
    pvect  presents momenta of particles: 
       4-momentum of  i^{th} particle ( i=0,1,...,nIn+nOut-1)  is 
       q[k]=pvect[4*i+k]  k=0,1,2,3; 
       q[0] - in particle energy, which always is positive.
       q[3] - specify projection of momentum on axis of collision.

    pName[i] (i=0,..nIn+nOut-1) contains name of i^th particle involved in
    reaction and pCode[i] is the corresponding PDG code. 

    Auxiliary functions which can help for construct usrFF are 
     
*/

extern int findval(char*name, double *value);
extern int qnumbers(char*pname, int *spin2, int * charge3, int * cdim);

/* The first one gives valueof  model parameter specified by its name.
If this  parameter indeed presented in the model then return value is zero 
and  parameter  'value' gets corresponding number. 

The qnumbers function gives  particle  quantum numbers: spin*2, 
(electric charge)*3 and dimension of color group representation. Return value
is PDG code which has to agree with pCode array data.
*/


/* Example: UMT(p1,p2) function which calculates transfer mass of 2 particles,
   for instance    UMT(e,Ne) - gives transverce mass of electron and neutrino.*/ 

double usrfun(char * name, int nIn, int nOut, double * pvect,char**pName,int*pCode)
{   char p1[10],p2[10];   // for 2 particles in MT(p1,p2)
    int i,j; 
    double sum=0;
    
    if(name==strstr(name,"MT("))   // name is started from "MT("
    {  //read p1&p2
       int np=sscanf(name+3,"%[^,],%[^)]",p1,p2);
       for(i=nIn;i<nIn+nOut;i++) 
       {  if(strcmp(p1,p2)==0) j=i+1;    /* if  p1==p2 */ else j=nIn; 
          for(  ;j<nIn+nOut;j++)
          if(strcmp(p1,pName[i])==0 &&  strcmp(p2,pName[j])==0) //find position of particles
          { double * q1=pvect+4*i, *q2=pvect+4*j;
            double Et1=sqrt(fabs(q1[0]*q1[0] - q1[3]*q1[3]) );    // transvers energy of the first particle
            double Et2=sqrt(fabs(q2[0]*q2[0] - q2[3]*q2[3]) );    // transvers energy of the second particle 
            sum+=sqrt( (Et1+Et2)*(Et1+Et2) -(q1[1]+q2[1])*(q1[1]+q2[1]) - (q1[2]+q2[2])*(q1[2]+q2[2])    ); // sqrt(E^2-PL^2)            
          }
       }      
    }  else { printf("Not defined user function %s\n",name); exit(2);}  
   
    return sum;
}
