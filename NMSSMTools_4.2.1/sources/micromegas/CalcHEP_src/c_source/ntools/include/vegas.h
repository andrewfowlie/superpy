#ifndef __VEGAS__
#define __VEGAS__


typedef  struct vegasGrid 
{
   
int ndim,     /* number of dimensions */ 
    ndmx;   
    long nCubes;
    double * x_grid;
    double * c_grid;
    float  * fMax;
} vegasGrid;
        
extern vegasGrid *  vegas_init
(int dim,  /* number of dimensions */
 int ndmx   /* size of grid */
);

extern void vegas_finish( vegasGrid * vegPtr);

extern int vegas_int(vegasGrid * vegPtr, 
 long ncall0,                       /* number of integrand calls */
 double alph,                       /* rate of grid improvement  */
 double(*fxn)(double *,double),     /* integrand */
 double *ti,                        /* integral estimation */ 
 double *tsi                        /* standard deviation */
);


extern int vegas_max(
vegasGrid * vegPtr, 
long  nCubes, 
long nRandom,
long nSimplex,
double milk,
double (*fxn)( double *,double,double*), 
double * eff
);


extern long vegas_events(
vegasGrid * vegPtr, 
long  nEvents,
double gmax, 
double (*fxn)( double *,double,double*), 
void (*out)(long ,int,char*,double*),
int recalc,   /* recalculate events in cube in case of new maximum */
double * eff,  /* efficiency */
int * nmax, /* max reached */
int * mult, /* partion of multiple events */
int * neg   /* partion of events with negative weght */
);

#endif
