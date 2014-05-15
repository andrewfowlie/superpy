/**************************************************************************
*    ccolor - c-program to colour weight calculation in SU(n=3) 
*    Program based on CompHEP version of the program.
*    Copyright (C)   1998, A.Kryukov (kryukov@theory.npi.msu.su)
*    URL:            http://theory.npi.msu.su/~kryukov
*--------------------------------------------------------------------------
*    Start           Jan 29, 1998
*    Version         1.2
*    Last revision   Feb. 06, 1998
*                    Mar. 12, 1998
*                    May  25, 1998  insert cmode arg. in calcCG
**************************************************************************/

#ifndef CCOLOR_H

#define CCOLOR_H

#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>

#include "physics.h"

/*################# definition from cweight.c #####################*/

#define NCOLOR   3    /*  SU(3)                           */
#define MAX_POW  32   /*  Maximum power of polinom of N_c */
#define MAXE     3    /*  Standard QCD                    */
#define MAXVT    3    /*                                  */
#define MAXGLEN  (2*MAXINOUT+2)    /*  16? Maximum length of CGraph        */
#define REVSP2T  2    /*  Sp(Ta*Tb)=1/RevSp2T*Delta(a,b)  */
#define CERRLEV1 0    /*  Run time error level            */
#define CERRLEV2 0    /*  Halt level                      */ 
#define CDEBLEV  00   /*  Total debug level               */ 


/*################################################################*/
   
typedef enum {zv=1,tv,g2,qg,g3} vtype;   /*  02/01/90                 */ 
                                      /* ZV    Zero vertex         */ 
                                      /* TV    Tranfer vertex      */ 
                                      /* G2    Transfer gluon vertex  */ 
                                      /* G3    Three gluon vertex  */ 
                                      /* QG    Qark-gluon vertex   */ 
                          /*  Range for C-graph length  */ 

typedef struct vertex 
   { 
      vtype        vt;           /* vertex type: zv, tv, g2, qg, g3 */
      int        e[MAXE];        /* array of edges (linked vrtex numbers) */
   }  vertex; 
   
typedef struct cgraph 
   {  long         n[MAX_POW];   /* Numerator of c-weight: polinom of N_c */ 
      long         d;            /* Denumerator           */ 
     long          r;            /* Power of N_c in denominator */
      int          en;           /* Name of next edge */
      int          gl;           /* Number of vertecies (graph length) */
      vertex       vl[MAXGLEN];  /* array of verticies */
   }  cgraph; 

typedef struct glist 
   {  cgraph       cg; 
      struct glist *next; 
   }  glist;

typedef struct weight 
   {  long        n[MAX_POW]; 
      long        d; 
      long        r;
      glist      *pgl; 
   }  weight; 


#include "cweight.h"

/*####################### end of defs ############################*/


#define CC_DEBUG 0

void initTokenCG(FILE *file);/* initiate reading process */
void initCG(weight *w);     /* initiate weight structure */
void readCG(weight *w);     /* read c-graph and build weight structure */
void printCG(cgraph *c);    /* print c-grapg */
void printPol(long *n,long d,long r);    /* print Polinom */
void calcCG(weight *w,int cmode);     /* calculate colour weight */

#endif /* CCOLOR_H */
