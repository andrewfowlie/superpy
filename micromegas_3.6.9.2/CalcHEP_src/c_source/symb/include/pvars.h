#ifndef __PVARS_
#define __PVARS_

#include"polynom.h"

typedef struct polyvars
   {
      int        nvar;   /* <============== Must be initialized to 0  */
      varinfo * vars;
   }  polyvars;

extern polyvars *vardef;

extern void  increaseVars(polyvars * v);
extern void  clearVars(polyvars * v);
extern void  addvar(char  * varname,  int    deg);
extern void  sortvar(void);
extern void  unite_vardef(polyvars *vardef_s,polyvars *vardef);
extern int   modelVarPos(char  * s);
extern int   scalarProductPos(int p1, int p2);

#define ALIG(n) (n>0 ?  8*((n-1)/8+1): 0 )

#define PPSHIFT 6
#endif
