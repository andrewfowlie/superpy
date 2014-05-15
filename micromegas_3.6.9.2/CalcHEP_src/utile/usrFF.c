#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<string.h>

extern double usrFF(int nIn, int nOut,  double * pvect,char**pName,int*pCode);
/*   
The usrFF function appears as factor at squared matrix element for 
Monte Carlo calculations.  CalcHEP code has a dummy version of this function 
which always return 1.  The dummy version is replaced on the  user one if 
its code is passed to CalcHEP linker via 'Libraries' model file. One can use 
CALCHEP and  WORK  environment variables  to specify path to the code. 
These variables are  automatically defined in calchep and calchep_batch  scripts. 
Also one can use any other environment variables defined separately. 
 Parameters of usrFF:
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


double usrFF(int nIn, int nOut, double * pvect,char**pName,int*pCode) 
{ retrun 1; }
