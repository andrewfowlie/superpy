#ifndef __ALPHAS2__
#define __ALPHAS2__
#include<stdio.h>
#include"nType.h"
extern int qcdmen_(void);
extern int w_alphaQCD(FILE *mode);
extern int r_alphaQCD(FILE *mode);
extern void i_alphaQCD(void);
extern int w_Scales(FILE *mode);
extern int r_Scales(FILE *mode);
extern void i_Scales(void);

extern void Scale(double*pv,double*qF,double*qR);
extern double (*sf_alpha)(double);
extern double  alpha_2(double Q);
#endif
