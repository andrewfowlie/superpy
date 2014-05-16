#include "micromegas.h"
#include "micromegas_aux.h"
#include "micromegas_f.h"

/*
numout*newProcess(char*Process)
{  int err;
   char lib[100];
   err=process2Lib(Process,lib);
   if(err) return NULL;      
   return getMEcode(0,ForceUG,Process,NULL,"",lib);
}
*/

void newprocess_(char*Process, int * address, int len1)
{ char cProcess[100];
  numout*cc;
  fName2c(Process,cProcess,len1);
  cc=newProcess(cProcess);
  if(!cc) *address=0;else memcpy(address,&cc,sizeof(cc)); 
}
  
     
void  procinfo1_(int*ccf, int *ntot, int * nin, int *nout)
{  numout*cc;
   memcpy(&cc,ccf,sizeof(cc));
   procInfo1(cc, ntot, nin, nout);}

void procinfo2_(int*ccf,int*nsub,char*name,double*mass,int len)
{ int ntot, nin, nout,i;
  numout*cc;
  char ** cname;  
  memcpy(&cc,ccf,sizeof(cc));
  procInfo1(cc, &ntot, &nin, &nout);
  cname=malloc((nin+nout)*sizeof(char*));
  
  procInfo2(cc,*nsub,cname, mass);
  for(i=0;i<nin+nout;i++) cName2f(cname[i],name+i*len,len);
}

