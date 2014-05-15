#include <stdio.h>
#include <stdio.h>
#include <math.h>   
#include <complex.h>

#include"nType.h"
        
extern int slhaReadW(char *fname,int mode);
extern int slhaRead(char *fname,int mode);
extern double slhaVal(char * Block, double Q, int nKey, ...);
extern int slhaValExists(char * Block, int nKey, ...);
extern double slhaWidth(int pNum);
extern int slhaWrite(char *fname);
extern int slhaWarnings(FILE*f);
extern int slhaDecayExists(int pNum);
extern double slhaBranch(int pNum,int N, int * nCh);
extern int initDiagonal(void);
extern int rDiagonal(int nDim,...); 
extern int rDiagonalA(int nDim,...);
extern REAL MassArray(int id,  int i);
extern REAL MixMatrix(int id, int i,int j); 
extern REAL MixMatrixU(int id, int i,int j);
extern int cDiagonalH(int Dim,...);
extern int cDiagonalA(int Dim,...);
extern int cDiagonalS(int Dim,...);
extern double complex cMixMatrix(int id,int i,int j);
extern double complex cMixMatrixU(int id,int i,int j);
extern int System(char * format, ...);
extern int openAppend(char * fileName);
extern int aPrintF(char * format,...);

extern double initQCD(double,double,double,double);
extern double initQCD5(double,double,double,double);
extern double MbEff(double);
extern double MtEff(double);
extern double McEff(double);
extern double alphaQCD(double);
extern double nfQCD(double Q);
extern double poleQmass(double M_M_, double alpha, int nf);

extern double complex HggF(double z);
extern double complex HggS(double z);
extern double complex HggV(double z);
extern double complex HggA(double z);

extern double complex Hgam1F(double z);
extern double complex Hgam1S(double z);
extern double complex Hgam1V(double z);
extern double complex Hgam1A(double z);
extern double Mbp(void);

#define  slhaVal0(block,scale)             slhaVal(block,scale,0)
#define  slhaVal1(block,scale,i1)          slhaVal(block,scale,1,i1)
#define  slhaVal2(block,scale,i1,i2)       slhaVal(block,scale,2,i1,i2)
#define  slhaVal3(block,scale,i1,i2,i3)    slhaVal(block,scale,3,i1,i2,i3)
#define  slhaVal4(block,scale,i1,i2,i3,i4) slhaVal(block,scale,4,i1,i2,i3,i4)
#define  slhaValExists0(block)             slhaValExists(block,0) 
#define  slhaValExists1(block,i1)          slhaValExists(block,1,i1)
#define  slhaValExists2(block,i1,i2)       slhaValExists(block,2,i1,i2)
#define  slhaValExists3(block,i1,i2,i3)    slhaValExists(block,3,i1,i2,i3)
#define  slhaValExists4(block,i1,i2,i3,i4) slhaValExists(block,4,i1,i2,i3,i4)

/* To avoid avto-prototyping  

extern double  sqrt(double);
extern double  sin(double);
extern double  cos(double);
extern double  tan(double);
extern double  asin(double);
extern double  acos(double);
extern double  atan(double);
extern double  exp(double);
extern double  log(double); 
extern double  pow(double,double);
extern double  fabs(double);
extern double  atan2(double,double);
extern double  log10(double);
extern double  sinh(double);
extern double  cosh(double);
extern double  tanh(double);
extern double  asinh(double);
extern double  acosh(double);
extern double  atanh(double);

extern double  creal(double complex);
extern double  cimag(double complex);
extern double  carg(double complex);
extern double  cabs(double complex);
extern double complex conj(double complex);
extern double complex cacos(double complex);
extern double complex casin(double complex);
extern double complex catan(double complex);
extern double complex ccos(double complex);
extern double complex csin(double complex);
extern double complex ctan(double complex);
extern double complex cacosh(double complex);
extern double complex casinh(double complex);
extern double complex catanh(double complex);
extern double complex ccosh(double complex);
extern double complex csinh(double complex);
extern double complex ctanh(double complex);
extern double complex cexp(double complex);
extern double complex clog(double complex);
extern double complex clog10(double complex);
extern double complex cpow(double complex,double complex);
extern double complex csqrt(double complex);
extern double complex cproj(double complexl);

extern int printf(char*, ...); 

extern double  if(double,double,double);

*/
