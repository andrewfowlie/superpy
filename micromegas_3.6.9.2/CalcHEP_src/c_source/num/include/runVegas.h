#ifndef __RUNVEGAS__
#define __RUNVEGAS__

#include<stdio.h>
#include"vegas.h"

extern int runVegas(void);
extern int runEvents(void);
extern void clearGrid(void);
extern void clearEventMax(void);
extern void newSession(void);
extern int nSess;

extern double inP1,inP2;

extern int saveVegasGrid( FILE * f); 
extern int readVegasGrid( FILE * f);

typedef struct vegas_integral
{  int itmx[2];
   long ncall[2];
   int freeze;
   double In, dI,khi2; 
   double s0,s1,s2;
   long  nCallTot;
   int n_it,old,tp;   
}vegas_integral;


extern vegas_integral integral; 

extern char * effInfo(void);

int saveEventSettings(FILE * f);
int readEventSettings(FILE * f);


#endif
