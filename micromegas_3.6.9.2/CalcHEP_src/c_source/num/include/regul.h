#ifndef __REGUL__
#define __REGUL__
#include "kinaux.h"
#include "chep_crt.h"


extern table regTab;

typedef struct {
    double pos;
    double width;
    int    power;
}  sing_struct;


extern double sngfun_(long *nn);
extern int wrtreg_(FILE *nchan);
extern int rdrreg_(FILE *nchan);
extern int getreg_(int *nsing, sing_struct *singar, double shift, double fmult,
int nfirst);
extern int fillRegArray(void);



typedef struct {
    char lvinvr[PLISTLEN];
    double rgmass, rgwdth;
    int nextrg, ndeg;
} invreg_;

extern invreg_ invreg_1[200];

#endif
