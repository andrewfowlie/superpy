#include<stdio.h>
#ifndef __4_VERTOR__
#define __4_VERTOR__

#include"nType.h"
extern REAL  vsqrt(REAL a);
extern REAL  vdot4(int i, int j,REAL*V);
extern void  vsum4(int i, int j, int k, int isg,REAL*V);
extern void  eps4(int n1, int n2, int n3, int n4,REAL*V);
extern void  pvFill(REAL  mass, REAL * mom, int pos,REAL*V); 
extern void  lvtonv(char * lv, int nin, int nvpos,REAL*V);
//extern REAL  pvect[400];

extern void incomkin(REAL m1, REAL m2, REAL p1, REAL p2,  
           REAL *sqrt_S_, REAL *Pcm_,REAL * rapidity_);
            
#endif
