#include"../../sources/micromegas.h"
#include"../../sources/micromegas_aux.h"
#include"pmodel.h"
#include"pmodel_f.h"

int readVarMSSM(char * fname)
{ int rdCode;
  char*vlist[35]={ "alfSMZ","MW","MZ","Ml","McMc","MbMb","Mtp","tb","MG1","MG2",
"MG3","Am","Al","At","Ab","Au","Ad","MH3","mu","Ml1","Ml2",
"Ml3","Mr1","Mr2","Mr3","Mq1","Mq2","Mq3","Mu1","Mu2","Mu3","Md1","Md2","Md3","wt"};

  rdCode = readVarSpecial(fname,35,vlist);

  return rdCode;
} 

int  readvarmssm_(char * f_name,int len)
{
  char c_name[100];
  fName2c(f_name,c_name,len);
  
  return readVarMSSM(c_name);
}
