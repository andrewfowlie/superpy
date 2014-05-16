#ifndef __CHESS_
#define __CHESS_

#include"sets.h"


typedef struct
{
   int weight, g5, vlnc;
   set ind;
   int link[6*maxvert];
}  vertinfostr ;

extern vertinfostr	vertinfo[6 * maxvert];
extern int  n_vrt;
extern int  prgcode[6 * maxvert][2];


extern void  makeprgcode(void);

#define MEMORY_OPTIM  0

#endif
