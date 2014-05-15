#include"../../sources/micromegas.h"
#include"pmodel.h"


void o1Contents(FILE * f)
{  
  fprintf(f,"\n~o1 = %.3f*bino %+.3f*wino %+.3f*higgsino1 %+.3f*higgsino2\n",
      findValW("Zn11"), findValW("Zn12"), findValW("Zn13"), findValW("Zn14"));
}
