#include "pmodel.h"
#include "../../sources/micromegas.h"
#include "../../sources/micromegas_aux.h"

/* ============================================================
   The function presented in this file are used only in calchep 
   interactive session. Thus the prototipes are never  checked.
   So we present them below only to avoid warnings at compilaton
   ============================================================= */
#undef deltarho 
#undef gmuon 
#undef bsgnlo 
#undef bsmumu 
#undef masslimits
                                                                                                 
extern forCalchep1 dMb,deltarho,gmuon,bsgnlo,bsmumu,masslimits,slha;
/*
extern double saveSM(double MbMb,double Mtp,double SW,double alfSMZ,
       double alfEMZ,double MZ,double Ml);
       
extern double  Omega(double qcdOk);
*/
/* =====  end of header part ========= */
/*
double saveSM(double MbMb,double Mtp,double SW,double alfSMZ,double alfEMZ,double MZ,double Ml)
{
  double Qmin;
  assignValW("MbMb",  MbMb); 
  assignValW("Mtp",   Mtp);
  assignValW("SW",    SW);
  assignValW("alfSMZ",alfSMZ);
  assignValW("alfEMZ",alfEMZ);
  assignValW("MZ",    MZ);
  assignValW("Ml",    Ml);
  Qmin=initQCD(alfSMZ,1.4,MbMb,Mtp);
  if(Qmin<0 || Qmin>MbMb) return 0; else return 1;
}
*/

double  dMb(double ok)     {return deltaMb();    }
double  deltarho(double ok){ return deltarho_(); } 
double  gmuon(double ok)   { return  gmuon_();   }
double  bsgnlo(double ok)  { return  bsgnlo_(NULL);  }
double  bsmumu(double ok)  { return  bsmumu_();  }
double  masslimits(double ok) {return masslimits_();}

/*
double  Omega(double qcdOk)
{ double omega;

  sortOddParticles(NULL);
  omega=darkOmega(NULL,1,1.E-6);

  return omega;
}
*/

extern int nSess;
double slha(double mssmOK)
{ char fname[20];
  sprintf(fname,"slha_%d.txt",nSess);
  slhaWrite(fname);
  return 1.;
}

extern double  saveSLHA(double x);
double  saveSLHA(double x){  return 1;}


static double SW,tb,MZ,mu;

double findGIparam(double ok)
{
  double Z[4][4], M[4],MM[4][4];
  int i,j,k,l;
  double mZero,mTot;
  double w,tw,cw,b,cb;
  double QSUSY=100;
 
  M[0]=slhaVal("MASS",QSUSY,1,1000022);
  M[1]=slhaVal("MASS",QSUSY,1,1000023);
  M[2]=slhaVal("MASS",QSUSY,1,1000025);
  M[3]=slhaVal("MASS",QSUSY,1,1000035);


 for(i=0;i<4;i++) for(j=0;j<4;j++) Z[j][i]=slhaVal("NMIX",QSUSY,2,i+1,j+1);


 mTot=0; 

 for( i=0;i<4;i++)
 { double s;

   for(l=0;l<4;l++)
   { for(s=0,k=0;k<4;k++) s+=Z[i][k]*M[k]*Z[l][k];
     MM[i][l]=s; 
     if(fabs(s)>mTot) mTot=fabs(s);
   }
 }

 mZero=0; 
 if(mZero<fabs(MM[0][1])) mZero=fabs(MM[0][1]);
 if(mZero<fabs(MM[2][2])) mZero=fabs(MM[2][2]);
 if(mZero<fabs(MM[3][3])) mZero=fabs(MM[3][3]);

 if(mZero>1.E-4*mTot) return 0;
 


 tw=-MM[0][2]/MM[1][2];
 w=atan(tw);
 SW=sin(w);
 cw=cos(w);
 tb=-MM[0][3]/MM[0][2];
 b=atan(tb);
 cb=cos(b);
 MZ= -MM[2][0]/cb/SW;
 mu=-MM[3][2]; 
 return 1;
}




double SWGI(double ok){return SW;}
double MZGI(double ok){return MZ;}
double tbGI(double ok){return tb;}
double muGI(double ok){return mu;}
