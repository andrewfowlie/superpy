#ifndef __READ_MDL_
#define __READ_MDL_

#include"file_scr.h"
extern table modelTab[5];

#define vars_tab modelTab[0]
#define func_tab modelTab[1]
#define prtcls_tab modelTab[2]
#define lgrng_tab modelTab[3]


extern int  readModelFiles(char * path, int l);
extern int  loadModel(int check, int ugForce);
extern int read2VarsParticles(void);
extern  int makeVandP(int rd ,char*path,int L, int mode,char*CaLCHEP);

#endif
