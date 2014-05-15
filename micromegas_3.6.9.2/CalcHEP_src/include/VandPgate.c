#include<string.h>
#include"VandP.h"


double usrFF_(int n_in, int n_out,double * pvect,char**pnames,int*pdg)
{ return  usrFF(n_in,n_out,pvect,pnames,pdg);}
double usrfun_(char * name,int n_in, int n_out,double * pvect,char**pnames,int*pdg) 
{ return  usrfun(name,n_in,n_out,pvect,pnames,pdg);}


extern double  alphaspdf_(double *Q );
extern void   getdatapath_(char* dirpath, int len);

extern void   initpdfsetbynamem_(int *P,char *name, int len);
extern void  evolvepdfm_(int* P,double *x,double *Q,double *f);
extern void  initpdfm_(int* P,int * nSet );
extern void  numberpdfm_(int* P,int * nMax);


double alpha_lha(double q ) { return alphaspdf_(&q); }

void   getdatapath(char* dirpath, int len){ getdatapath_(dirpath, len);}
void   initpdfsetbynamem(int *P,char *name, int len){ initpdfsetbynamem_(P,name,len);}
void  evolvepdfm(int* P,double *x,double *Q,double *f){evolvepdfm_(P,x,Q,f) ;}
void  initpdfm(int* P,int * nSet,double*xMin,double*xMax,double*qMin,double*qMax)
{
  extern void getxmaxm_(int*,int*,double *);
  extern void getxminm_(int*,int*,double *);
  extern void getq2maxm_(int*,int*,double *);                     
  extern void getq2minm_(int*,int*,double *);

  initpdfm_(P,nSet ) ;
  getxmaxm_(P,nSet,xMax); 
  getxminm_(P,nSet,xMin);
  getq2maxm_(P,nSet,qMax);  *qMax=sqrt(fabs(*qMax));
  getq2minm_(P,nSet,qMin);  *qMin=sqrt(fabs(*qMin));
//  printf("limits: %E %E %E %E\n", *xMin,*xMax,*qMin,*qMax);
}

void  numberpdfm(int* P,int * nMax){ numberpdfm_(P,nMax);}

extern int findval(char *name,double *val);
extern int qnumbers(char*pname, int *spin2, int * charge3, int * cdim);

int findval(char *name,double *val)
{ int i;
  for(i=0;i<nModelVars+nModelFunc;i++) 
  if(!strcmp(name,varNames[i])){ *val= varValues[i]; return 0;}  
  return 1;
}


int qnumbers(char*pname, int *spin2, int * charge3, int * cdim)
{
   int n,sign;
   for(n=0;n<nModelParticles;n++)
   { 
     if(!strcmp(pname,ModelPrtcls[n].name )) {sign=1; break;} 
     if(!strcmp(pname,ModelPrtcls[n].aname)) {sign=-1;break;}
   }
   if(n==nModelParticles) return 0;

   if(spin2)   *spin2  =ModelPrtcls[n].spin2;
   if(charge3) *charge3=sign*ModelPrtcls[n].q3;
   if(cdim)    
   {  *cdim   =ModelPrtcls[n].cdim; 
      if(sign==-1 &&(*cdim==3 || *cdim==-3)) (*cdim)*=-1;
   }
   return sign*ModelPrtcls[n-1].NPDG;
}
