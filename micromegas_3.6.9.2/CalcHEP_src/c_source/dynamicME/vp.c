#include<stdio.h>
#include<string.h>
#include<math.h>

#include"vp.h"
#include"dynamic_cs.h"

int pTabPos(char * name)
{
  int i;
  for(i=0;i<nModelParticles;i++)
  { 
    if(!strcmp(name,ModelPrtcls[i].name )) return   i+1;
    if(!strcmp(name,ModelPrtcls[i].aname)) return -(i+1);
  }
  return 0;
}

char * pdg2name(int pdg)
{
  int i;
  if(pdg==0) return NULL;

  for(i=0;i<nModelParticles;i++)
  {          if(ModelPrtcls[i].NPDG==pdg) return ModelPrtcls[i].name;
     else  { if(ModelPrtcls[i].NPDG==-pdg) return ModelPrtcls[i].aname;}
  }   
  {
    static char name[20];
    sprintf(name,"#%d\n",pdg);
    return name;
  }
  return NULL;
} 

double pMass(char * name)
{
  char *nm;
  int n=pTabPos(name);
  if(!n){printf("Wrong particle name '%s'\n",name); return 0;}
  nm=ModelPrtcls[abs(n)-1].mass;
  if(nm[0]=='0') return 0; else 
  { REAL *ma=varAddress(nm);
    return fabs(*ma);
  }
}

int pNum(char * name)
{
  int n=pTabPos(name);
  if(!n){printf("Wrong particle name ''%s''\n",name); return 0;}
  if(n>0)  return  ModelPrtcls[abs(n)-1].NPDG;
  else     return -ModelPrtcls[abs(n)-1].NPDG;
}

int qNumbers(char*pname, int *spin2, int * charge3, int * cdim)
{
  int n=pTabPos(pname);
  int sign=1;
  int pdg;
  if(!n) { if(pname[0]=='#' && sscanf(pname+1,"%d",&pdg)==1) return pdg; else return 0;}
  if(n<0){ n=-n; sign=-1;}
  if(spin2)   *spin2  =ModelPrtcls[n-1].spin2;
  if(charge3) *charge3=sign*ModelPrtcls[n-1].q3;
  pdg=sign*ModelPrtcls[n-1].NPDG;
  if(cdim)    
  { *cdim   =ModelPrtcls[n-1].cdim; 
    if(sign==-1 &&(*cdim==3 || *cdim==-3)) (*cdim)*=-1;
  }
  return pdg;
}

REAL * varAddress(char *name)
{int i;
 for(i=0;i<nModelVars+nModelFunc;i++) if(!strcmp(name,varNames[i]))return varValues+i;
 return NULL;
}

int  findVal(char * name, double * val)
{
  int i;
  for(i=0;i<nModelVars+nModelFunc;i++)
  { 
    if(strcmp(name,varNames[i])) continue;
    *val=varValues[i] ;
    return 0;
  }
  return 2;
}

double  findValW(char*name)
{ double val;
  if(findVal(name,&val)) {printf(" %s not found\n",  name); return 0;}
  else return val;
}
  
