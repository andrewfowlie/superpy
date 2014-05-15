#ifndef __DIAGRAMS__
#define __DIAGRAMS__

#include"model.h"

/*========== Pukhov representation ================== */

#define ldiagram (2 * MAXINOUT - 3)
/*  maximum number of particles in diagram  */

typedef   particleNumType  decayDiagram[ldiagram]; 

typedef  struct   adiagram
{
    decayDiagram    dgrm0;
    char             delMark;
    int nsub;
} adiagram;

typedef struct csdiagram
   {
      decayDiagram  dgrm1,dgrm2;
      int        lnk[MAXINOUT];
      int           mult;
      unsigned          del; 
      char          status;   /* -2-outOfMemory,-1-deleted,
                                 0-Rest, 1-calculated,2-Zero  */
      int nsub;                           
   }  csdiagram;

#define maxvert (MAXINOUT - 2)
/*  maximal # of verteces in amplitude           */
#define nullvert 253    
/*# of next vert for  unused edge              */

/*================== Taranov representation =========== */

#define IN_PRTCL   1
#define OUT_PRTCL  2
#define PLR_PRTCL  4

typedef struct vertlink
{  
  int  vno, edno; /* # of vert, # of edge in vert (slot) */
}  vertlink;

typedef struct edgeinvert
   {
      int          lorentz;
      int          moment;
      int          prop;
      int          partcl;
      vertlink     nextvert;
      vertlink     link;
   }  edgeinvert;

typedef edgeinvert vert0[MAXVALENCE];


typedef struct vampl 
{ 
   int          size, outno;       /*  how many  verts and outgoing edges  */ 
   int          valence[maxvert];
   vertlink     outer[MAXINOUT];   /*  adresses of external edges */ 
   vert0        vertlist[maxvert]; /*  array of verts  */ 
} vampl; 
            

typedef struct vcsect
   {
      int         sizel, sizet;
      long        symnum, symdenum, clrnum, clrdenum;
      int          /* 1..4 */ valence[2 * maxvert];
      vert0        vertlist[2 * maxvert];
   }  vcsect;


extern vcsect     vcs;

extern void  transfdiagr(csdiagram  * diag,  vcsect *     vcs);
extern void  mkverts(decayDiagram diag1,vampl* vlist1);

extern void InOutPrtclsNumb( decayDiagram a, int * numb, int sort);
extern void proccessName(decayDiagram a, char * txt );

extern void  decompose( vcsect vcs,  vampl *left,  vampl * right);

extern void printDiagram(vampl * vlist);
extern void printCsDiagram(vcsect * vlist);
#endif
