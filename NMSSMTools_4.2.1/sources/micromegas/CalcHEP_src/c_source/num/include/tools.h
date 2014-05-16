#ifndef __TOOLS__
#define __TOOLS__

extern double dinter_(double x, int n, double *xi, double *yi);
extern double gammai_(int n, double a);
extern double convol_(double(*f1)(double),double(*f2)(double), double b1, double b2, double x, double eps);
extern double divy_(double );

#endif
