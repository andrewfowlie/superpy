#include"../../sources/micromegas.h"
#include"../../sources/micromegas_aux.h"
#include"pmodel.h"
#include"pmodel_f.h"

int readVarRHNM(char * fname)
{
  char*vlist[44]={"EE","alfSMZ","SW","s12","s23","s13","Mm","Ml","Mu","Md","McMc","Ms",
  "MbMb","Mtpole","MZ","MH","wtp","wZ","wW","v","rc","g10","gH","gZ","gZp","gUp","gtl",
  "gll","gtr","gbr","grl","xWr","gHZ","MLZP","MZp","MWp","cr","Mul","Mur","Mn5","Mtl","Mbl","Mbr","mixzzp"};

  return readVarSpecial(fname,44,vlist);
} 

int  readvarrhnm_(char * f_name,int len)
{
  char c_name[100];
  fName2c(f_name,c_name,len);
  
  return readVarRHNM(c_name);
}
