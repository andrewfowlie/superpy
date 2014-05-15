#include"../../sources/micromegas_aux.h"
#include"pmodel.h"


void o1Contents(FILE * f)
{ double ok=0; 
  fprintf(f,"\n~o1 = ");

  fprintf(f,"(%.3f%+.3f*i)*bino+",     findValW("Zn11r"),findValW("Zn11i")); 
  fprintf(f,"(%.3f%+.3f*i)*wino+",     findValW("Zn12r"),findValW("Zn12i")); 
  fprintf(f,"(%.3f%+.3f*i)*higgsino1+",findValW("Zn13r"),findValW("Zn13i")); 
  fprintf(f,"(%.3f%+.3f*i)*higgsino2", findValW("Zn14r"),findValW("Zn14i")); 
  fprintf(f,"\n");
}

int readVarCPVMSSM(char * fname)
{
  char*vlist[43]={
  "alfSMZ","Mtp", "MbMb","McMc","EE",  "SW",  "Ml", "MHc", "aMu","fiMu",
  "aM1",   "aM2", "aM3", "fiM1","fiM2","fiM3","Ml2","Ml3", "Mr2","Mr3",
  "aAt",   "fiAt","aAb", "fiAb","aAl", "fiAl","aAm","fiAm","aAu","fiAu",
  "aAd",   "fiAd","Mq2", "Mq3", "Mu2", "Mu3", "Md2","Md3", "wt", "MZ",
  "tb",    "aAe", "fiAe"};

  return readVarSpecial(fname,43,vlist);
} 

