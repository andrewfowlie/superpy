#ifndef __SAVERES_
#define __SAVERES_

#include"physics.h"
#include"polynom.h"
#include"denominators.h"


extern  denom_struct   denom[2 * maxvert - 2];

extern int  denrno;

extern void  saveanaliticresult(poly rnum, poly factn, poly factd, vcsect vcs, int nFile);

#endif
