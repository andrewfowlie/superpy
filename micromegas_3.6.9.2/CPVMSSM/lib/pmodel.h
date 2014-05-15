#ifndef __CPSH__
#define __CPSH__

#ifdef __cplusplus
extern "C" {
#endif

#include<stdio.h>

#include"../../CalcHEP_src/c_source/SLHAplus/include/SLHAplus.h"


extern void o1Contents(FILE * f);
 

extern double cpHiggs(double,double,double,double,double,double,double,double,
 double,double,double,double,double,double,double,double,double,double,double,
 double,double,double,double,double,double,double,double,double,double,double,
 double,double,double,double,double,double);

extern int readVarCPVMSSM(char * fname);
extern int loopGamma(double * cs_gz, double *cs_gg);

#ifdef __cplusplus
}
#endif 

#endif
