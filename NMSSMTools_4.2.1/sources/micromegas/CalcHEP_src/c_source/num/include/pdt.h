#ifndef __PDT__
#define __PDT__

#include <stdio.h>
#include <string.h>
#include<stdlib.h>
#include<math.h>

typedef struct pdtList
{ struct pdtList * next; 
  char *  name;      /* title name of distribution             */     
  long  beamP;     /* MC-number of beam particles            */
  char * file;       /* name of file where it is stored        */
  int  * partons;    /* ordering number of parton in the list of functions*/
  int  * items;
} pdtList; 

typedef struct pdtStr
{ double mass;       /* mass of composite particle             */
  long beamP;
  long parton;
  int  nq;           /* number of points in Q-scale grid       */ 
  int  nx;           /* number of points in X-scale grid       */
  double * x_grid;   /* data for the X-grid                    */
  double * q_grid;   /* log() data for the Q-grid              */
  double * x_grid_aux;
  double * q_grid_aux;
  
  double * alpha;    /* data for QCD-alpha(Q) corresponding to */ 
                     /* the Q-grid points*/ 
  double * strfun;   /* data for interpolation of distibution corresponding */
  double * strfun_aux;
  double * aux;
  double(*interpolation)(double, double,struct pdtStr *);
                     /* to  (Q-grid)*(X-grid) points;          */ 
                     /* index of (X-grid) is rinning first     */
  double * q_grid_cteq;
  double   Q0_cteq;
  double x_min;
  double q_min;
  double q_max;      /* boundaries for x and q arguments       */
  double pow0;       /* factor x^pow *(1-x)^pow1  must be      */   
  double pow1;       /* applied after interpolation            */ 
  double q_threshold;/* threshold q for heavy quarks           */
  long nSmallX;
  long nSmallQ;
  long nLargeX;
  long nLargeQ;
  
  int    qt0;        /* position of first q point above threshold*/
  int  approx;
} pdtStr;  


extern long   makePdtList(char * file, pdtList ** list);
                     /* scan 'file' to detect does it contain  */
                     /* data for 'parton'. If yes, the         */
                     /* corresponding information is added to 'list' */

extern void   delPdtList(pdtList * list);
                     /* free  memory allocated for 'list' */
                     
extern int    getPdtData(char * file, int n_parton, pdtStr * data );
                     /* read 'file' and  fill items of  'data' */ 

extern void   freePdtData( pdtStr * data);
                    /* free  memory allocated for data items   */
                    /* and assigne NULL to  them.              */

extern double interFunc(double x, double q, pdtStr * W); 
                     /* interpolates data for given x and q    */ 
                     /* according to information stored in W.  */
                     /* result should be multiplied by the     */  
                     /* power factors later on                 */
extern double interAlpha(double q, pdtStr * W );
                     /* interpolates data for QCD-alpha(Q)     */
                     
extern int checkPartons( int * pNum, pdtList * L);                     

#endif
