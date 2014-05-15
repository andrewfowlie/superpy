#include"pmodel.h"
#include"../../sources/micromegas_aux.h"
#include"pmodel_f.h"
#include"pmodel.h"

int nmhwarn_(int *file)
{ 
  char fname[100];
  FILE*f;
  int r;

  if(*file==0) return slhaWarnings(NULL);
                                                                                
  sprintf(fname,"%d.tmptxt",getpid());
  f=fopen(fname,"w");
  r=slhaWarnings(f);
  fclose(f);
  fortreread_(file,fname,strlen(fname));
  unlink(fname);
  return r;
}


void o1contents_(int *Nch)
{
  char fname[20];
  FILE*f;

  sprintf(fname,"%d.tmptxt",getpid());
  f=fopen(fname,"w");
  o1Contents(f);
  fclose(f);
  fortreread_(Nch,fname,strlen(fname));
  unlink(fname);
}

int nmssmewsb_(int *mode ){ return nmssmEWSB(*mode );}

int  nmssmsugra_(double *m0, double* mhf, double* a0, double* tb,
double*sgn, double*Lambda, double *aLambda, double*aKappa)
{
  return  nmssmSUGRA(*m0, *mhf, *a0, *tb, *sgn, *Lambda,*aLambda, *aKappa);
} 

int leshinput_(char * fname, int len)
{ char c_name[200];
  fName2c(fname, c_name,len);
  return  lesHinput(c_name);
}

int  readvarnmssm_(char * f_name,int len)
{
  char c_name[100];
  fName2c(f_name,c_name,len);
  
  return readVarNMSSM(c_name);
}

int hbblocks_(char*fname, int len)
{  char cname[100];
   fName2c(fname,cname,len);
   return   HBblocks(cname);
}  
 