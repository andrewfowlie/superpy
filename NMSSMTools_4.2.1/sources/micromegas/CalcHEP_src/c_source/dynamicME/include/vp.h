#ifndef __VP__
#define __VP__

#include "../../../include/VandP.h"

extern int      pTabPos(char * name);
extern double   pMass(char * name);
extern int      pNum(char * name);
extern int      qNumbers(char *pname, int *spin2, int * charge3, int * cdim);
extern char *   pdg2name(int pdg);
extern REAL*    varAddress(char * name);
 

#endif
