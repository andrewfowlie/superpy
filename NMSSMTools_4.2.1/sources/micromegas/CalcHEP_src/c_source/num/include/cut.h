#ifndef __CUT__
#define __CUT__

#include<stdio.h>
#include "file_scr.h"
#include "kinaux.h"
#include "phys_val.h"

extern int fillCutArray(void);
extern int rdrcut_(FILE *);
extern int wrtcut_(FILE *);

extern double  calcCutFactor(double*V);

extern int rancor_(double *vmin, double *vmax, double shift, double fmult, int n);

extern void  rancor_t(double * cosmax, double hsum , double fmult, double Ecm,
       double mq, double pcmtilda,  double mp1, double ptMin);

extern table cutTab;


typedef struct {
       int aux;
       char title[50];
       char key[4];
       physValRec *pLists;
       int minon, maxon;
       double  cvmin, cvmax;
} invcut_;


 
extern invcut_ invcut_1[60];
extern int nCuts;

#endif
