#ifndef __NUM_IN_
#define __NUM_IN_

#include<stdlib.h>
#include<string.h> 
#include<math.h>
#include "nType.h"

extern double alpha_2(double);
typedef REAL (DNN)(double, REAL *,int *);
typedef REAL (FNN)(double,REAL*,REAL*,COMPLEX*,REAL*);

extern  double Fmax;
extern  REAL Helicity[2];
extern  REAL HelicityN[2];
extern  int    CalcConst;
extern  int    indx_(int k,int l);
extern  void   sprod_(int ntot, REAL * momenta, REAL*DP);
extern  int    prepDen(int nden, int nin, double BWrange2,   
                       REAL * dmass,REAL * dwidth, char **q,REAL * mom, 
                       REAL *Q0,COMPLEX*Q1,REAL*Q2);
#endif
