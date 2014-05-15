#include"../../sources/micromegas_aux.h"
#include"pmodel_f.h"
#include"pmodel.h"


void o1contents_(int *file)
{
  char fname[20];
  FILE*f;

  sprintf(fname,"%d.tmptxt",getpid());
  f=fopen(fname,"w");
  o1Contents(f);
  fclose(f);
  fortreread_(file,fname,strlen(fname));
  unlink(fname);
}

int  readvarcpvmssm_(char * f_name,int len)
{
  char c_name[100];
  fName2c(f_name,c_name,len);
  
  return readVarCPVMSSM(c_name);
}

int loopgamma_(double * cs1, double *cs2){ return loopGamma(cs1,cs2);}

