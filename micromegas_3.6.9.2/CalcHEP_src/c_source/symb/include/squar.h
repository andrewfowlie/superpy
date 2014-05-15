#ifndef __SQUAR_
#define __SQUAR_

#include"model.h"

typedef int permut[MAXINOUT];

typedef struct permListStr
{
  struct permListStr * next;
  permut      perm;
} permListStr;

typedef struct permListStr *permlist;

typedef struct amplListStr
   {
      struct amplListStr * next;
      decayDiagram      dgrm;
      permut       perm;
      permlist     gen;
      unsigned         dim;
   }  amplListStr;
typedef struct amplListStr *ampllist;



extern int squaring(void);

extern ampllist  labeledDiagrams(FILE * diagrp);
extern void clearDiagrams(ampllist a);


#endif
