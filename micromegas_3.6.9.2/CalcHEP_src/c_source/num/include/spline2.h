#ifndef __SPLINE2__
#define __SPLINE2__

extern void progonca(int N,double *f,double *d);
extern void b_and_c(int N,double *f,double *b,double *c,double *d);
extern double spline_for_graph(double iks,double *f,double *b,double *c,double *d);
extern void SPLINE2(int koeff,int N,double *Iexp,double *dIexp,double *I);

#endif
