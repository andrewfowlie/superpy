#ifndef __PROCESS_
#define __PROCESS_

#include"syst2.h"
#include"model.h"
#include"physics.h"

#define whohowMAX 300
typedef struct whohow
{  int    who, how;
}  whohow[whohowMAX];

typedef struct hadron
{
   char        name[P_NAME_SIZE];
   shortstr    contents;
   int         pow;
   int         parton[1000];
   int         polarized[1000];
}  hadron;

extern whohow     liminsp, LimQ;
extern whohow     limout;
extern void  nilprtcl(whohow      p_list);


extern shortstr processch;
char limpch[STRSIZ], deloutch[STRSIZ];
extern int  nin, nout, n_x;   /* Number of X-particles */
extern int  enter(void);
extern  hadron hadrons[MAXINOUT];

extern int polarized(int p, int Prtcl);

extern int ZWmax,ZWmin;

#endif
