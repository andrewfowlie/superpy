/*
 Copyright (C) 1997, Slava Ilyin 
*/

#include"sf_lsr.h"
#include<math.h>
#include<string.h>

int p_lsr__(int * pNum) { if(pNum[0]==22 && pNum[1]==0) return 1; else return 0;}

int mc_lsr__(int i) {return 22;}
void n_lsr__(int i, char *name) { strcpy(name,"Laser photons");}

int r_lsr__(int i, char *name)
{  if (strcmp(name, "Laser photons") != 0) return 0; else return 1; }

int i_lsr__(int i, double* be, double * mass) { *be=1.; *mass=0.; return 1;}

int m_lsr__(int i,int*pString) {return 1;}

double c_lsr__(int i, double x, double q)
{

    static int first = 1;
    double x0 = 4.82;
    double xmax=x0/(1+x0);
    static double rnorma;
 
    if (first) 
    {  double d__2= x0 + 1.;
        first = 0;
	rnorma = (1. - 4. / x0 - 8. / (x0 * x0)) * log(x0 + 1.) + .5 + 8. 
		/ x0 - 1. / (d__2 * d__2 * 2.);
    }
    if (x > xmax)   return  0.; else 
	return (1. - x + 1. / (1. - x) * (1. - x * 4. / x0 * (1. - x / (
		x0 * (1. - x))))) / rnorma;
}

