/*-------------------------------------------------*
 * C -version by Victor Edneral, Moscow University *
 *        The latest revision of 23.01.1994        *
 *-------------------------------------------------*/

#ifndef __PLOT__
#define __PLOT__

extern void   plot_1(double xMin, double xMax, int dim, 
                    double *f, double *ff,char* upstr, char* xstr, char* ystr);
extern void plot_2(double hMin1,double hMax1,int nBin1,
                   double hMin2,double hMax2,int nBin2,
            double * f,double *df,char *proces,char* xname,char * yname);
#endif
