#ifndef __SPLINE__
#define __SPLINE__

extern  void spline3make(int n, double *f, double *c);

extern double spline3val(int n, double * f, double * c, double x, double * deriv);

extern void spline3smooth(double k, double p, int N, double *Iexp,double *dIexp, double *I);

#endif
