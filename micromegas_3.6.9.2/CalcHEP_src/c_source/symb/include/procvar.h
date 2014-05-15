#ifndef __PROCVAR_
#define __PROCVAR_

#include "physics.h"

#define VAR_NAME_SIZE_EXT  20 

typedef struct
   {
      char    alias[VAR_NAME_SIZE_EXT];
      double  tmpvalue;
      int     num; 
      int     used;
   }  singlevardescription;

extern singlevardescription *vararr;
extern int nProcessVar;

extern int  initvararray(int nsub, char key,int width);

#endif
