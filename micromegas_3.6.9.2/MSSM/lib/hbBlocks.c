#include"../../sources/micromegas.h"
#include"../../sources/micromegas_aux.h"
#include "pmodel.h"

#define SQR(x) (x)*(x)
int  HBblocks(char * fname)
{ FILE * f=fopen(fname,"a");
  double tb,sb,cb,alpha,sa,ca,ta,samb,camb,dMb,MbHl,MbSM,MbH,MbH3,Q;
  if(!f) return 1;
  Q=findValW("Q");
  if(slhaDecayExists(pNum("h")) <0)  slhaDecayPrint("h", f);
  if(slhaDecayExists(pNum("H")) <0)  slhaDecayPrint("H", f);
  if(slhaDecayExists(pNum("H3"))<0)  slhaDecayPrint("H3",f);
  if(slhaDecayExists(pNum("t")) <0)  slhaDecayPrint("t", f);
  if(slhaDecayExists(pNum("H+"))<0)  slhaDecayPrint("H+",f);

  tb=findValW("tB");  
    sb=tb/sqrt(1+tb*tb);
    cb=1/sqrt(1+tb*tb);
  alpha=findValW("alpha");
    sa=sin(alpha);
    ca=cos(alpha);
    ta=sa/ca;
    samb=sa*cb-ca*sb;
    camb=ca*cb+sa*sb;
  dMb=findValW("dMb");


  MbSM=findValW("Mb");
  MbH= MbSM/(1+dMb)*(1+dMb*ta/tb);  
  MbH3=MbSM/(1+dMb)*(1-dMb/tb/tb);
  MbHl=MbSM/(1+dMb)*(1-dMb/ta/tb);
 
  fprintf(f,"Block HiggsBoundsInputHiggsCouplingsBosons\n");
  fprintf(f,"# Effective coupling normalised to SM one and squared\n");
  fprintf(f,"# For (*) normalized on Sin(2*W)\n"); 
  fprintf(f," %12.4E  3    25    24    24 # higgs-W-W \n",        SQR(samb)  );
  fprintf(f," %12.4E  3    25    23    23 # higgs-Z-Z \n",        SQR(samb)  );
  fprintf(f," %12.4E  3    25    25    23 # higgs-higgs-Z \n",    0.   );

  { assignVal("Q",pMass("h"));
    calcMainFunc();
    fprintf(f," %12.4E  3    25    21    21 # higgs-gluon-gluon\n",  SQR(findValW("LGGh")/findValW("LGGSM")) );           
    fprintf(f," %12.4E  3    25    22    22 # higgs-gamma-gamma\n",  SQR(findValW("LAAh")/findValW("LAASM")) );
  }                          
 
  fprintf(f," %12.4E  3    35    24    24 # higgs-W-W \n",        SQR(camb)  );
  fprintf(f," %12.4E  3    35    23    23 # higgs-Z-Z \n",        SQR(camb)  );
  fprintf(f," %12.4E  3    35    25    23 # higgs-higgs-Z \n",    0.  );
  fprintf(f," %12.4E  3    35    35    23 # higgs-higgs-Z \n",    0.  );
  
  
  { assignVal("Q",pMass("H"));
    calcMainFunc();
    fprintf(f," %12.4E  3    35    21    21 # higgs-gluon-gluon\n",SQR(findValW("LGGH")/findValW("LGGSM"))  );   
    fprintf(f," %12.4E  3    35    22    22 # higgs-gamma-gamma\n",SQR(findValW("LAAH")/findValW("LAASM"))  );
  }
  

  fprintf(f," %12.4E  3    36    24    24 # higgs-W-W \n",        0.  );
  fprintf(f," %12.4E  3    36    23    23 # higgs-Z-Z \n",        0.  );
  
  { assignVal("Q",pMass("H3"));
    calcMainFunc();
    fprintf(f," %12.4E  3    36    21    21 # higgs-gluon-gluon\n",SQR(findValW("LGGH3")/2/findValW("LGGSM")) );
    fprintf(f," %12.4E  3    36    22    22 # higgs-gamma-gamma\n",SQR(findValW("LAAH3")/2/findValW("LAASM")) );             
  }
  
  fprintf(f," %12.4E  3    36    25    23 #*higgs-higgs-Z \n",    SQR(camb)  );
  fprintf(f," %12.4E  3    36    35    23 #*higgs-higgs-Z \n",    SQR(samb)  );
  fprintf(f," %12.4E  3    36    36    23 #* higgs-higgs-Z \n",   0.  );

  fprintf(f,"Block HiggsBoundsInputHiggsCouplingsFermions\n");
  fprintf(f,"# Effective coupling normalised to SM one and squared\n");
  fprintf(f," %12.4E   %12.4E   3    25     5    5 # higgs-b-b \n"    ,SQR((sa/cb)*(MbHl/MbSM)),0.);
  fprintf(f," %12.4E   %12.4E   3    25     6    6 # higgs-top-top \n",SQR(ca/sb)              ,0.);
  fprintf(f," %12.4E   %12.4E   3    25    15   15 # higgs-tau-tau \n",SQR(sa/cb)              ,0.);

  fprintf(f," %12.4E   %12.4E   3    35     5    5 # higgs-b-b \n"    ,SQR((ca/cb)*(MbH/MbSM))  ,0.);
  fprintf(f," %12.4E   %12.4E   3    35     6    6 # higgs-top-top \n",SQR(sa/sb)              ,0.);  
  fprintf(f," %12.4E   %12.4E   3    35    15   15 # higgs-tau-tau \n",SQR(ca/cb)  ,0.);

  fprintf(f," %12.4E   %12.4E   3    36     5    5 # higgs-b-b \n"    ,0.,SQR(tb*(MbH3/MbSM)));
  fprintf(f," %12.4E   %12.4E   3    36     6    6 # higgs-top-top \n",0.,SQR(1/tb)          );
  fprintf(f," %12.4E   %12.4E   3    36    15   15 # higgs-tau-tau \n",0.,SQR(tb)            );
  
  assignValW("Q",Q);
  calcMainFunc();    
  fclose(f);
   
  return 0;
}
