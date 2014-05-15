#ifndef  __SIMPSON__
#define __SIMPSON__
extern double gauss( double (*func)(double),double a,double b, int n);
extern double gauss345(double (*func)(double), double a, double b, double eps,int * err_code);
extern double simpson( double (*func)(double),double a,double b, double  eps);
#endif
