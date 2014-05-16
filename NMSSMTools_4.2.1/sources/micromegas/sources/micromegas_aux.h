#ifndef __MICRO_AUX__
#define __MICRO_AUX__

#include"micromegas.h"
#include"../CalcHEP_src/include/VandP.h"

#ifdef __cplusplus
extern "C" {
#endif 

#include "../CalcHEP_src/c_source/dynamicME/include/dynamic_cs.h"
#include "../CalcHEP_src/c_source/dynamicME/include/vp.h"


/*==================== displayPlot =================*/
extern void displayPlot(char*title,char*xName,char*yName,double xMin,double xMax,int dim,double*f,double*ff);

/*======================
   Lock/UnLock service
========================*/
extern int checkLockFile(int *delay);
extern void removeLock(int fd);
/*
extern char *  prepareWorkPlace(void);
extern int cleanworkPlace(void);
*/
extern char*WORK;

/*=======================
      Odd particles
========================*/
extern int Nodd;
extern  ModelPrtclsStr * OddPrtcls;
extern int createTableOddPrtcls(void);
extern char * OddParticles(void);
extern char * EvenParticles(void);     
/*=============================================
C->F and F->C    string convertation     
===============================================*/

extern void  cName2f(char*c_name,char*f_name,int len);
extern void  fName2c(char*f_name,char*c_name,int len);

/*=============================================
Fortran output
===============================================*/
extern void fortreread_(int* N, char * fname, int len);

/*====================================
   Tools for integration 
=====================================*/

extern double simpson(double (*func)(double), double a, double b, double eps);
extern double gauss(double (*func)(double), double a, double b, int n);
extern int  odeint(double*ystart,int nvar,double x1,double x2, double eps,
                   double h1, void(*derivs)(double,double*,double *));


/*==== Tool  for interpolation  ====*/
extern double polint2(double x, int n,  double *xa, double *ya);
extern double  polint3(double x, int n,  double *xa, double *ya);
extern double  polint4(double x, int n,  double *xa, double *ya);
extern int buildInterpolation( double (*Fun)(double), double x1,double x2, double eps,
                                int * N, double ** xa, double **ya);
/*======= special functions ========*/
extern double bessk0(double x);
extern double bessk1(double x);
extern double bessk2(double x);
extern double K2pol(double x); /*bessk1(1/x)*exp(1/x)*sqrt(2/M_PI/x);*/
extern double K1pol(double x); /*bessk2(1/x)*exp(1/x)*sqrt(2/M_PI/x);*/

/*======== read Table =============*/

extern int readTable(char * fileName, int *Ncolumn, double **tab);


/* Hidden interface with CalcHEP */
extern int FError;
extern int OnlyTEQ0;

extern int pname2lib(char*pname, char * libname);
extern numout*newProcess_(int twidth, int model,int UG,char*Process,
          char * excludeVirtual,char*excludeOut,char*lib,int usr);


//=====  2->2 processes ===========

extern double (*sqme22)(int nsub, double GG, double *pvect, int * err_code); 

extern int     kin22(double PcmIn,double * pmass);
extern double  dSigma_dCos(double  cos_f);
extern int  nsub22;

#define NTOF(X) extern forCalchep1 X; double X(double ok){return findValW(#X);}

typedef  double (forCalchep1)(double);
typedef  double (forCalchep2)(double,double);

/*  Loop integrals I1 ... I5 */

extern double   LintIk(int i,double MSQ,double MQ,double MNE);

extern int readVarSpecial(char *fname, int nVar, char ** names);

extern double parton_x( int pNum, double  Q);
extern double parton_alpha(double q);
extern double parton_distr(int pNum, double x, double q);



extern double convStrFun2(double x, double q, int pc1, int pc2, int ppFlag ); 
                           /* result of convolution of structure functions 
                              of pc1 and pc2 particles  */

extern int  wimpPos(void);


extern int vPolar(int out1,int out2,int out3,double*left,double*right,double*lng);
extern double  plazmaWidth(char *process,double T);

extern double cs23MM(numout*cc, int nsub, double Pcm, int fast,int ii3, double M45min, double M45max);


extern double amoeba(double *p, double * y, int ndim, double (*f)(double *), double eps, int *nCalls);
                                                    


#include"../CalcHEP_src/include/num_in.h"

#include"microPath.h"


#ifdef __cplusplus
}
#endif 


#endif
