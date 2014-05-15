#ifndef __IS_ISR__
#define __IS_ISR__

#include<stdio.h>

extern int  p_isr__(int *pNum);
extern void n_isr__(int i, char *name);
extern int  r_isr__(int i, char *name);
extern int  m_isr__(int i,int*);
extern int  mc_isr__(int i);
extern int i_isr__(int i,double* be, double * mass);
extern double c_isr__(int i, double x,double q);

#endif
