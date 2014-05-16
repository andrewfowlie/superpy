#ifndef __COMPOSITES__
#define __COMPOSITES__

#include "file_scr.h"

extern int fillCompositeArray(void);
extern int rdrcomp_(FILE *);
extern int wrtcomp_(FILE *);


extern table compTab;
extern int nComps,nCompParts[60];
extern char compName[60][4], compParts[60][60][4];


#endif
