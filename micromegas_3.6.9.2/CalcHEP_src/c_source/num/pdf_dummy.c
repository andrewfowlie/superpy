#include <string.h>
#include"lha.h"


double  alphaspdf_(double *Q ){ return 0.12;}

void   getdatapath_(char* dirpath, int len)
{ int i;
  for(i=0;i<len;i++) dirpath[i]=' ';
}

void   initpdfsetbynamem_(int *P,char *name, int len)
{ return; }

void  numberpdfm_(int* P,int * nMax) { *nMax=0;}
void  evolvepdfm_(int* P,double *x,double *Q,double *f)
{  int i; for(i=0;i<13;i++) f[i]=0;}
void  initpdfm_(int* P,int * nSet ){ return; }

void getxmaxm_(int*P,int*N,double *xMax){ *xMax=1;}
void getxminm_(int*P,int*N,double *xMin){ *xMin=0;}
void getq2maxm_(int*P,int*N,double *qMax){ *qMax=1;}                  
void getq2minm_(int*P,int*N,double *qMin){ *qMin=1.E10;}
