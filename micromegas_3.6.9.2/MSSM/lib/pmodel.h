#ifndef  __MSSM__
#define  __MSSM__

#ifdef __cplusplus
extern "C" {
#endif 

#include"../../CalcHEP_src/c_source/SLHAplus/include/SLHAplus.h"
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<unistd.h>


/*=====================================================
   MSSM Parameters input at low scale
=======================================================*/

extern int  suspectEwsbMSSM(void);
extern int  isajetEwsbMSSM(void);

/*=============================================
  MSSM Parameters motivated by SUGRA scenario
  =============================================*/
extern int  suspectSUGRA(
 double tb,  double gMG1,double gMG2,double gMG3,
 double gAl, double gAt, double gAb, double sgn, double gMHu, double gMHd,
 double gMl1,double gMl3,double gMr1,double gMr3,
 double gMq1,double gMq3,double gMu1,double gMu3,double gMd1,double gMd3
                        );
extern int isajetSUGRA(
 double tb, double gMG1,double gMG2,double gMG3,
 double gAl, double gAt, double gAb, double sgn, double gMHu, double gMHd,
 double gMl1,double gMl3,double gMr1,double gMr3,
 double gMq1,double gMq3,double gMu1,double gMu3,double gMd1,double gMd3
                      );


extern int  softSusySUGRA(
 double tb, double gMG1,double gMG2,double gMG3,
 double gAl, double gAt, double gAb, double sgn, double gMHu, double gMHd,
 double gMl1,double gMl3,double gMr1,double gMr3,
 double gMq1,double gMq3,double gMu1,double gMu3,double gMd1,double gMd3
                         );

extern int sphenoSUGRA(
 double tb, double gMG1,double gMG2,double gMG3,
 double gAl,double gAt, double gAb, double sgn, double gMHu, double gMHd,
 double gMl1,double gMl3,double gMr1,double gMr3,
 double gMq1,double gMq3,double gMu1,double gMu3,double gMd1,double gMd3
                      );
/*=============================================  
  SUGRA with fixed mu and MH3  
 =============================================*/
 extern int  suspectSUGRAnuh( double tb, double gMG1,double gMG2,double gMG3,double gAl, double gAt, double gAb,
 double gMl1,double gMl3,double gMr1,double gMr3,double gMq1,double gMq3,double gMu1,double gMu3,double gMd1,double gMd3,
 double mu,double MA );

 extern int isajetSUGRAnuh( double tb, double gMG1,double gMG2,double gMG3,double gAl, double gAt, double gAb,
 double gMl1,double gMl3,double gMr1,double gMr3,double gMq1,double gMq3,double gMu1,double gMu3,double gMd1,double gMd3,
 double mu,double MA );

 extern int  softSusySUGRAnuh( double tb, double gMG1,double gMG2,double gMG3,double gAl, double gAt, double gAb, 
 double gMl1,double gMl3,double gMr1,double gMr3,double gMq1,double gMq3,double gMu1,double gMu3,double gMd1,double gMd3,
 double mu,double MA );

 extern int sphenoSUGRAnuh(
 double tb, double gMG1,double gMG2,double gMG3, double gAl,double gAt, double gAb,
 double gMl1,double gMl3,double gMr1,double gMr3,double gMq1,double gMq3,double gMu1,double gMu3,double gMd1,double gMd3,
 double mu,double MA );

/*=============================================
  MSSM Parameters motivated by other  scenarios
  =============================================*/

extern int suspectAMSB(double am0,double m32,double stgbeta,double ssgnmu0);
extern int isajetAMSB(double am0,double m32,double tb,double sgn);
extern int softSusyAMSB(double am0,double m32,double tb,double sgn);
extern int sphenoAMSB(double am0,double m32,double tb,double sgn);

extern int isajetGMSB(double Lambda, double Mmess,double tb, double SGNMU,double N5, double cGrav,
double Rsl, double dmH_d2, double dmH_u2, double d_Y, double n5_1, double n5_2, double n5_3);



/*===============================================
     Les Houches interface
 ================================================*/
extern int  readLesH(char*fname, int SM );
extern int  writeLesH(char*fname);     

int lesHinput(char * fname);



/*================================================
  Experimental constrains  on the  MSSM parameters  
  ===============================================*/
extern double deltarho_(void);
extern double gmuon_(void);
extern double gmuonold_(void);
extern int    masslimits_(void);
extern double bsgnlo_(double *SMbsg);
extern double bsmumu_(void);
extern double btaunu_(void);
extern double deltaMb(void);
extern double deltaMs(void);
extern double deltaMd(void);
extern double deltaMl(void);
extern double deltaMc(void);
extern double Rl23_(void);
extern double dtaunu_(double *dmunu);
extern int dMQcorrections;
extern int loopGamma(double * cs_gz, double *cs_gg);
extern int callSuperIsoSLHA(void);

#define deltarho   deltarho_
#define gmuon      gmuon_
#define masslimits masslimits_
#define bsgnlo     bsgnlo_
#define bsmumu     bsmumu_
#define btaunu     btaunu_
#define Rl23       Rl23_
#define dtaunu     dtaunu_

/*======================
  printMasses + 
=========================*/

extern void o1Contents(FILE * f);

/*====readVarMSSM====*/

extern int readVarMSSM(char * fname);

/*===== calculation of box diagrams =====*/
extern double qbox_(void);
      
extern int MSSMDDtest(int loop, double*pA0,double*pA5,double*nA0,double*nA5);


extern int  HBblocks(char * fname); 


#ifdef __cplusplus
}
#endif 

#endif
