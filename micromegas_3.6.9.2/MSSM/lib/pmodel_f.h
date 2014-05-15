#ifndef  __MSSM_F_
#define  __MSSM_F_


/*=====================================================
   MSSM Parameters input at low scale
=======================================================*/

extern int  suspectewsbmssm_( int * lCon);
/* 
  integer function suspectEwsbMSSM(lcOn)
  integer lsOn
*/

extern int  isajetewsbmssm_(int *lCon);
/*
  integer function isajetEwsbMSSM(lCon)
  integer lCon
*/

/*=============================================
  MSSM Parameters motivated by SUGRA scenario
  =============================================*/
extern int  suspectsugra_(
double*tb,double*MG1,double*MG2,double*MG3,double*Al,double*At,double*Ab,
double*sgn,double*MHu,double*MHd,double*Ml1,double*Ml3,double*Mr1,double*Mr3,
double*Mq1,double*Mq3,double*Mu1,double*Mu3,double*Md1,double*Md3 );
/*
     integer function  suspectSUGRA(tb,MG1,MG2,MG3,Al,At,Ab,sgn,MHu,MHd,
   > Ml1,Ml3,Mr1,Mr3,gMq1,Mq3,Mu1,Mu3,Md1,Md3)
     real*8 tb,MG1,MG2,MG3,Al,At,Ab,sgn,MHu,MHd,Ml1,Ml3,Mr1,Mr3,gMq1,
   >  Mq3,Mu1,Mu3,Md1,Md3                           
*/

extern int  isajetsugra_(
double*tb,double*MG1,double*MG2,double*MG3,double*Al,double*At,double*Ab,
double*sgn,double*MHu,double*MHd,double*Ml1,double*Ml3,double*Mr1,double*Mr3,
double*Mq1,double*Mq3,double*Mu1,double*Mu3,double*Md1,double*Md3 );
/*
     integer function  isajetSUGRA(tb,MG1,MG2,MG3,Al,At,Ab,sgn,MHu,MHd,
   > Ml1,Ml3,Mr1,Mr3,gMq1,Mq3,Mu1,Mu3,Md1,Md3)
     real*8 tb,MG1,MG2,MG3,Al,At,Ab,sgn,MHu,MHd,Ml1,Ml3,Mr1,Mr3,gMq1,
   >  Mq3,Mu1,Mu3,Md1,Md3                           
*/

extern int  softsusysugra_(
double*tb,double*MG1,double*MG2,double*MG3,double*Al,double*At,double*Ab,
double*sgn,double*MHu,double*MHd,double*Ml1,double*Ml3,double*Mr1,double*Mr3,
double*Mq1,double*Mq3,double*Mu1,double*Mu3,double*Md1,double*Md3 );
/*
     integer function  softSusySUGRA(tb,MG1,MG2,MG3,Al,At,Ab,sgn,MHu,MHd,
   > Ml1,Ml3,Mr1,Mr3,gMq1,Mq3,Mu1,Mu3,Md1,Md3)
     real*8 tb,MG1,MG2,MG3,Al,At,Ab,sgn,MHu,MHd,Ml1,Ml3,Mr1,Mr3,gMq1,
   >  Mq3,Mu1,Mu3,Md1,Md3                           
*/

extern int  sphenosugra_(
double*tb,double*MG1,double*MG2,double*MG3,double*Al,double*At,double*Ab,
double*sgn,double*MHu,double*MHd,double*Ml1,double*Ml3,double*Mr1,double*Mr3,
double*Mq1,double*Mq3,double*Mu1,double*Mu3,double*Md1,double*Md3 );
/*
     integer function  sphenoSUGRA(tb,MG1,MG2,MG3,Al,At,Ab,sgn,MHu,MHd,
   > Ml1,Ml3,Mr1,Mr3,gMq1,Mq3,Mu1,Mu3,Md1,Md3)
     real*8 tb,MG1,MG2,MG3,Al,At,Ab,sgn,MHu,MHd,Ml1,Ml3,Mr1,Mr3,gMq1,
   >  Mq3,Mu1,Mu3,Md1,Md3                           
*/

/*=============================================
  MSSM Parameters motivated by other  scenarios
  =============================================*/
extern int suspectamsb_(double*am0,double*m32,double*stgbeta,double*sgn);
/*
    integer function suspectAMSB(am0,m32,tg,sgn)
    real*8 am0,m32,tg,sgn
*/

extern int isagetamsb_(double*am0,double*m32,double*stgbeta,double*sgn);
/*
    integer function isajetAMSB(am0,m32,tg,sgn)
    real*8 am0,m32,tg,sgn
*/

extern int softsusyamsb_(double*am0,double*m32,double*stgbeta,double*sgn);
/*
    integer function softsusyAMSB(am0,m32,tg,sgn)
    real*8 am0,m32,tg,sgn
*/

extern int sphenoamsb_(double*am0,double*m32,double*stgbeta,double*sgn);
/*
    integer function sphenoAMSB(am0,m32,tg,sgn)
    real*8 am0,m32,tg,sgn
*/

/*===============================================
      Les Houches accord interface
 =================================================*/
extern int readlesh_(char *fname , int*SM,  int len );
/* 
     integer function readLesH(fname,SM)
     character*(*) fname
     int SM
*/  

extern int writelesh_(char *fname, int len );
/*
     integer function  writeLesH(fname)
     character*(*) fname
*/
               

   

/*================================================
  Experimental constrains  on the  MSSM parameters  

 real*8  function  deltarho()
 real*8  function  gmuon()
 integer function  masslimits()
 real*8  function  bsgnlo()
 real*8  function  bsmumu()

 ===============================================*/

void o1contents_(int *N);
/*
  subroutine  o1Contents(int *N )
*/



extern double deltamb_(void);
/* 
    real*8 function deltaMb()
*/

extern int readvarmssm_(char * name, int len);
/*
     integer function readVarMSSM(name)
     character*(*)name     
*/       
#endif

