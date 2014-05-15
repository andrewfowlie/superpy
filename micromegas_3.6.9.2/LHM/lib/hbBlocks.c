#include"../../sources/micromegas.h"
#include"../../sources/micromegas_aux.h"
#include "pmodel.h"


#define SQR(x) (x)*(x)
int  HBblocks(char * fname)
{ FILE * f=fopen(fname,"w");
  double Q;
  if(!f) return 1;
  Q=findValW("Q");
  
 fprintf(f,"Block Mass\n 25  %E # Higgs Mass\n\n",findValW("MH"));
   
  slhaDecayPrint("H",f);
  slhaDecayPrint("t",f);
  slhaDecayPrint("~H+",f);

  fprintf(f,"Block HiggsBoundsInputHiggsCouplingsBosons\n");
  fprintf(f,"# Effective coupling normalised to SM one and squared\n");
  fprintf(f,"# For (*) normalized on Sin(2*W)\n"); 
  fprintf(f," %12.4E  3    25    24    24 # higgs-W-W \n",       SQR(( 1.+findValW("del")*findValW("del")/12.)*
( 1.-findValW("del")*findValW("del")*findValW("vh")*findValW("vh")/findValW("v")/findValW("v")/3.)) );
  fprintf(f," %12.4E  3    25    23    23 # higgs-Z-Z \n",       SQR(( 1.+findValW("del")*findValW("del")/12.)*
( 1.-findValW("del")*findValW("del")*findValW("vh")*findValW("vh")/findValW("v")/findValW("v")/3.)) );
  fprintf(f," %12.4E  3    25    25    23 # higgs-higgs-Z \n",   0. );
 
  { assignVal("Q",pMass("H"));
    calcMainFunc();
    fprintf(f," %12.4E  3    25    21    21 # higgs-gluon-gluon\n",  SQR(findValW("LGGH")/findValW("LGGSM")) );           
    fprintf(f," %12.4E  3    25    22    22 # higgs-gamma-gamma\n",  SQR(findValW("LAAH")/findValW("LAASM")) );
  }      


  fprintf(f,"Block HiggsBoundsInputHiggsCouplingsFermions\n");
  fprintf(f,"# Effective coupling normalised to SM one and squared\n");
  fprintf(f," %12.4E   %12.4E   3    25     5    5 # higgs-b-b \n"    ,SQR(-1+findValW("del")*findValW("del")*findValW("vh")*findValW("vh")/findValW("v")/findValW("v")/4.) ,0.);
  fprintf(f," %12.4E   %12.4E   3    25     6    6 # higgs-top-top \n",SQR(1+findValW("mtcorr")/4/findValW("v")/findValW("v")*findValW("B00014")),0.);
  fprintf(f," %12.4E   %12.4E   3    25    15   15 # higgs-tau-tau \n",1.,0.);

  fclose(f);

  assignValW("Q",Q);
  calcMainFunc();  
   
  return 0;
}

int  hbblocks_(char * fname,int len)
{
  char * cname=malloc(len+2);
  int err;
  fName2c(fname,cname,len);
  err=HBblocks(cname);
  free(cname);
  return err;
}
            