#ifndef __PREPDIAG_
#define __PREPDIAG_

#include"sets.h"
#include"ghosts.h"

typedef char momsum[MAXINOUT+1];

typedef struct
{
  int  len;
  int  g5;
  int  vv[2*maxvert], ll[2*maxvert], nint[2*maxvert],intln[2*maxvert][2];
/*  char invrt[2*maxvert]; */
  int  spin[2*maxvert];
  int  ind[2*maxvert][2];
} fermloopstp;


typedef	 struct
{
  algvertptr   lgrnptr;
  arr4byte     subst;
  int          r_vert;
}  vertexhlp;

typedef  struct{int  vrt1, ln1, vrt2, ln2;} linkhlp;

extern vertexhlp    vertexes[2*maxvert];
extern linkhlp      massindpos[5*maxvert];
extern fermloopstp  fermloops[maxvert];
extern set    setmassindex,setmassindex0;
extern int    nloop;
extern int    consLow;
extern int    fermmap[2*maxvert];
extern char   inoutmasses[MAXINOUT][VAR_NAME_SIZE];
extern momsum momdep[3*maxvert];

extern void   preperdiagram(void);
extern void   coloringvcs(hlpcsptr  currentghst);
extern void   attachvertexes(void);
extern void   findReversVert(void);
extern set   findmassindex(void);
#endif
