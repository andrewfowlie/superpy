#include<math.h>
#include"pmodel.h"
#include"pmodel_aux.h"
#include"pmodel_f.h"
#include"../../sources/micromegas.h"



int readLesH(char *fname, int mode)
{ /* mode 0 - EWSB;  1 - SUGRA/AMBS; 2 - FILE */
  int err;  
  err=slhaRead(fname,0);
  FillVal(2);
  if(err) 
  {  printf("Problem with file format. Probably '%s' is not SLHA format file\n",fname);
     
  }   
  if(slhaWarnings(NULL)) 
  {  printf("slhaRead  warnings:\n");
     slhaWarnings(stdout);
  }   
  if(err) return err;

  return 0;
}

int lesHinput(char * fname)
{ int err= readLesH(fname, 2);
  if(err==0) FillVal(2);
  return err;
}
