#ifndef __ALPHA__
#define __ALPHA__

#include<stdio.h>

extern double alpha(int nf, int odr, double lambda,  double dscale);

extern double findLambda(int nf,int odr, double alpha0, double M);

extern int writeAlpha(FILE*f,int nf,int ordr,double lambda, int nfMx,
                      double Mc,double Mb,double Mt,int N, double *q);

#endif 
