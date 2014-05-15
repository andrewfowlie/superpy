#ifndef __PARAM__
#define __PARAM__

#include"nType.h"

extern int checkParam(void);

extern int selectParam(int x,int y,  char*mess, void ** pscrPrt,
   int for22, int rd_on, int vars_on,  int func_on, REAL ** varPos, char * varName,int*nPos);

extern int change_parameter(int x,int y, int for22);

extern void show_depend(int x, int y);

extern REAL Pcm22;

#endif
