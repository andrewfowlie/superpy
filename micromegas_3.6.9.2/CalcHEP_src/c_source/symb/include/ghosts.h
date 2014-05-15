#ifndef __GHOSTS_
#define __GHOSTS_

#include"model.h"
#include"diagrams.h"

#define vbosonmark 0
#define ghostmark 1
#define antighostmark 2
#define sbosonmark 3
#define xbosonmark (-1)
#define ybosonmark (-2)


typedef char  hlpcsect[2 * maxvert][MAXVALENCE];
#define hlpcsptr struct hlpcsrec *
typedef struct hlpcsrec
   {
      hlpcsect     hlpcs;
      char     sgn;
      unsigned         num, maxnum;
      hlpcsptr     next;
   }  hlpcsrec;
#undef hlpcsptr
typedef struct hlpcsrec *hlpcsptr;

#define j_lgrptr struct j_lgrrec *
typedef struct j_lgrrec
   {
      j_lgrptr     next;
      algvertptr   lgrnptr;
   } j_lgrrec;
#undef j_lgrptr
typedef struct j_lgrrec *j_lgrptr;

extern j_lgrptr j_lgrarray[2 * maxvert];

extern void  generateghosts(vcsect   * vcs, hlpcsptr * alll);
extern void  ghostsForAmplitude(vampl * ampl, hlpcsptr * ghosts);
extern void  eraseghosts(hlpcsptr  gst);
extern void  vertinlgr(arr4byte     vert, int  nvert, arr4byte subst,
					                     algvertptr * lgr);

#endif
