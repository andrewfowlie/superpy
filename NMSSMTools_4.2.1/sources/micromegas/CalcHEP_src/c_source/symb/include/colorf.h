#ifndef __COLORS__
#define __COLORS__

#include"diagrams.h"

#define MAXE     3                 /*  Standard QCD                    */


typedef enum {zv=1,tv,g2,qg,g3} vtype;/*  02/01/90                 */
                                      /* ZV    Zero vertex         */
                                      /* TV    Tranfer vertex      */
                                      /* G2    Transfer gluon vertex  */
                                      /* G3    Three gluon vertex  */
                                      /* QG    Qark-gluon vertex   */
                          /*  Range for C-graph length  */

typedef struct cvertex
   {
      vtype        vt;           /* vertex type: zv, tv, g2, qg, g3 */
      int        e[MAXE];        /* array of edges (linked vrtex numbers) */
   }  cvertex;


typedef struct factor
   { 
      long * nc;
      int  len;
      int  dpow;
      long dc;  /* Denumerator           */
   }  factor;


extern factor * colorFactor (int nv, cvertex * vl);
extern void fct_num_calc(factor * fct,int Nc, long * n, long *d);
extern void fct_print(factor * fct, char *s);

extern void t2k2(vcsect* g, int * nv, cvertex * vl);
extern vtype typev(vert0 v,int valence);
#endif  
