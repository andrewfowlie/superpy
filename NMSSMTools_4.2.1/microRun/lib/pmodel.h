#ifndef __NMSSM__
#define __NMSSM__

#include<stdio.h>

#ifdef __cplusplus
extern "C" {
#endif


extern void   o1Contents(FILE * f);
extern int    nmssmEWSB(void); 
extern int    nmssmSUGRA(double m0, double mhf,   double a0, double tb,
                         double sgn,double Lambda,double aLambda, double aKappa);
extern int    readVarNMSSM(char *fname);
extern int    readSLHA(char * fname);
extern int    NMHwarn(FILE * f);

extern double  bsgnlo_(double *M, double*P);
extern double  deltamd_(double *M, double*P);
extern double  deltams_(double *M, double*P);
extern double  bsmumu_(double *M, double*P);
extern double  btaunu_(double *M, double*P);
extern double  gmuon_(double *M, double*P);
extern int     HBblocks(char * fname);
extern int     loopGamma(double * csAA, double *csAZ);

#define bsgnlo     bsgnlo_
#define deltaMd    deltamd_
#define deltaMs    deltams_
#define bsmumu     bsmumu_
#define btaunu     btaunu_
#define gmuon      gmuon_

#ifdef __cplusplus
}
#endif 

#endif
