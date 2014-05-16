#ifndef __NUM_OUT_int
#define __NUM_OUT_int

#include<stdlib.h>
#include<string.h> 
#include<math.h>

#include"nType.h"

extern  int nin_int;
extern  int nout_int;
extern  int nprc_int;
extern  int nvar_int;
extern  int nfunc_int;

extern char * (*pinf_int)(int nsub, int nprtcl,  REAL* pmass, int*pnum);
extern int (*pinfAux_int)(int nsub, int nprtcl, int* spin2, int*color,int*neutral);
extern char ** polarized_int;
extern char ** varName_int;

extern double (*sqme_int)(int nsub,double GG, REAL * momenta, int * err);
extern int (*calcFunc_int)(void);
extern int *twidth_int, *gswidth_int, *gtwidth_int;
extern double *BWrange_int;
extern REAL *va_int;
extern char*hiddenf;

extern void (*build_cb_int)(int nsub); 
extern void (*destroy_cb_int)(void);    
extern int *cb_pow_int;   
extern int *cb_nc_int; 
extern int ** cb_chains_int;
extern REAL ** cb_coeff_int;


#define  DENOMINATOR_ERROR   2
#endif

