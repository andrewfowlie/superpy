#include<stdio.h>
#include<stdlib.h>
#include"micromegas.h"
#include"micromegas_aux.h"

//#define PRINT

#define Gconst (0.7426E-30) /* m/g */
#define Vlight (299792458.) /* m/s */  
#define mp_g   (1.673E-24)  /* proton mass in grams */
#define mp_gev (0.939)   
#define GeVfm  (5.067731162)

/*====================  CAPTURE RATES ===============================*/


static double Rcm, *rTab,*tTab,*rhoTab, *phiTab, *aFraction;
static int nTab;


static double RS=6.96E10;   /* cm */

static double (*fvStat)(double)=NULL;
static double v_stat,MA, muX, FFalpha;

/*
double Tcapture(double T, double w, double v)
{ double mu=Mcdm/MA;
  double muP= (mu+1)/2;
  double muM= (mu-1)/2;
  double aP=sqrt(MA/(2*T))*(muP*v + muM*w);
  double aM=sqrt(MA/(2*T))*(muP*v - muM*w);  
  double bP=sqrt(MA/(2*T))*(muM*v + muP*w);
  double bM=sqrt(MA/(2*T))*(muM*v - muP*w);  
   

 return - T/(MA*mu*mu)*(  mu*(-aP*exp(-aM*aM) - aM*exp(-aP*aP))/sqrt(M_PI) 
                         +(mu/2 -mu*aP*aM - muP*muM)*( erf(aP)-erf(-aM)) + muP*muP*( erf(bP)-erf(-bM))*exp(-Mcdm/(2*T)*(v*v-w*w)));  

}
*/

static numout* cc23; 
extern double cs23(numout*cc, int nsub, double Pcm, int ii3);


static double r3Integrand(double r3)
{
  double  r,vescQ,rhoGeVcm3,mfrac,w2rate;
  r=pow(r3,0.333333333333333);
  vescQ=2*polint2(r/Rcm,nTab,rTab,phiTab);
  w2rate= vescQ-v_stat*v_stat*muX;
  if(w2rate<=0) return 0;
  rhoGeVcm3=polint2(r/Rcm,nTab,rTab,rhoTab)/mp_g;   // in proton mass units
  mfrac=polint2(r/Rcm,nTab,rTab,aFraction);
  if(FFalpha==0) return rhoGeVcm3*mfrac*w2rate; 
  else  return  rhoGeVcm3*mfrac*(1+muX)*( exp(-v_stat*v_stat*FFalpha) - exp(-(v_stat*v_stat+vescQ)*FFalpha/(1+muX)))/FFalpha;
}

static double vIntegrand(double v)
{ double u=fvStat(v); 
  double phiv,Rmax;
  if(!u)return 0;

  v_stat=v*1.E3/Vlight;  // km -> c units 
  phiv=muX*v_stat*v_stat/2;

  if(phiv> phiTab[0]) return 0;
  if(phiv <= phiTab[nTab-1]) Rmax=Rcm; else   Rmax=Rcm*polint2(phiv,nTab,phiTab,rTab); 
  return u*(4*M_PI/3.)*simpson(r3Integrand,0,Rmax*Rmax*Rmax,1.E-4)*Vlight*Vlight*1.E-1; // in cm
}



static double sIntegrand(double r3)
{
  double  r,rhoGeVcm3,mfrac;
  r=pow(r3,0.333333333333333);
  rhoGeVcm3=polint2(r/Rcm,nTab,rTab,rhoTab)/mp_g;   // in proton mass units
  mfrac=polint2(r/Rcm,nTab,rTab,aFraction);
  return rhoGeVcm3*mfrac;   
}

static void derivePhi(double r, double *x, double *dx)
{  double r2=r*r;
   dx[0]=4*M_PI*r2*polint2(r/Rcm,nTab,rTab,rhoTab);
   dx[1]=Gconst*100*x[0]/r2;
//printf(" r=%E dx[0]=%e dx[1]=%E\n", r,dx[0],dx[1]);    
}


static void fillGravPotential(void)
{ double x[2],phi0,r1,r0,m1;
  int i;
  r0=rTab[0]*Rcm;
  if(r0==0)
  {  phiTab[0]=0;
     r1=rTab[1]*Rcm;
     x[0]=4./3.*M_PI*rhoTab[0]*r1*r1*r1;
     phiTab[1]=x[1]=2./3.*M_PI*rhoTab[0]*r1*r1*Gconst*100;
     i=2;
  } else 
  { 
    x[0]=4./3.*M_PI*rhoTab[0]*r0*r0*r0;        
    phiTab[0]=x[1]=2./3.*M_PI*rhoTab[0]*r0*r0*Gconst*100;
//printf("phiTab[0]=%E\n",phiTab[0]);    
    i=1;
  }  
  for( ;i<nTab;i++)
  {
    odeint(x,2,Rcm*rTab[i-1], Rcm*rTab[i],1.E-3,Rcm*(rTab[i]-rTab[i-1])/3,derivePhi);
    phiTab[i]=x[1];       
  }    
  phi0=x[0]*Gconst*100/Rcm;
//printf("VescSurface=%E\n", sqrt(2*phi0)*Vlight); 
  for(i=0;i<nTab;i++) phiTab[i]= phiTab[nTab-1]-phiTab[i]+phi0;
//printf("VescCenter=%E\n", sqrt(2*phiTab[0])*Vlight);
}



#define SUNPOINTS 1268 
static double  rSun[SUNPOINTS], tSun[SUNPOINTS], rhoSun[SUNPOINTS], Hsun[SUNPOINTS],He4sun[SUNPOINTS],
He3sun[SUNPOINTS],C12sun[SUNPOINTS],N14sun[SUNPOINTS],O16sun[SUNPOINTS],phiSun[SUNPOINTS];


#define EARTHPOINTS 101

static double OEarth[EARTHPOINTS],NaEarth[EARTHPOINTS],MgEarth[EARTHPOINTS],AlEarth[EARTHPOINTS],
SiEarth[EARTHPOINTS],PEarth[EARTHPOINTS],SEarth[EARTHPOINTS],CaEarth[EARTHPOINTS],CrEarth[EARTHPOINTS],
FeEarth[EARTHPOINTS],NiEarth[EARTHPOINTS];

static double rEarth[EARTHPOINTS],tEarth[EARTHPOINTS],rhoEarth[EARTHPOINTS],phiEarth[EARTHPOINTS];


static int readSunData(void)
{
  static int rdOK=0;
  int i;
  if(!rdOK) 
  { FILE*f;
    char fname[300];
    sprintf(fname,"%s/sources/data_nu/SSM.dat",micrO);
    
    f=fopen(fname,"r"); 
    fscanf(f,"%*[^\n]"); 
    for(i=0;i<SUNPOINTS;i++)
    {
      fscanf(f,"%*lf %lf %lf %lf %*lf %*lf %lf %lf %lf %lf %lf %lf ",
      rSun+i, tSun+i, rhoSun+i, Hsun+i,He4sun+i, He3sun+i, C12sun+i,N14sun+i,O16sun+i);
//    printf("%E %E %E %E %E %E %E %E %E\n", rsun[i], tsun[i], rhosun[i], Hsun[i],He4sun[i], He3sun[i], C12sun[i],N14sun[i],O16sun[i]);
//      printf("frac= %E\n",(1-( Hsun[i]+He4sun[i]+ He3sun[i]+ C12sun[i]+N14sun[i]+O16sun[i]))/O16sun[i] ); 
    }
    fclose(f);
  }    
  rTab=rSun;
  tTab=tSun;
  rhoTab=rhoSun;
  phiTab=phiSun;
  nTab=SUNPOINTS;
  Rcm=RS; 

  if(!rdOK) {fillGravPotential(); rdOK=1;}
  return 0;
}

static double captureSun(double(*vfv)(double),  double pA0, double nA0, double pA5, double nA5)
{ 
  /*                 H    He      C     N     O    Ne     Mg    Si    S     Fe    Na    Al   Cl    Ar     Ca   Cr   Ni */
  int      A[17] ={  1,   4   , 12   ,14   ,   16,   20,   24,   28,   32,   56,   23,   27,   35,   40,   40,  52,  59  };
  int      Z[17] ={  1,   2   ,  6   , 7   ,    8,   10,   12,   14,   16,   26,   11,   13,   17,   18,   20,  24,  28  };
  double P10[17] ={ 12,  10.93,  8.39, 7.78, 8.66, 7.84, 7.53, 7.51, 7.14, 7.45, 6.17, 6.37, 5.50, 6.18, 6.31, 5.64, 6.23};  
//  double P10[17] ={ 12,  10.93,  8.39, 7.78, 8.66, 8.08, 7.58, 7.55, 7.33, 7.50 , 6.33, 6.47, 5.50, 6.40, 6.36, 5.64, 6.25};  
  double  sumI,csSDp,mu;
  double * sunFractions[5]={Hsun,He4sun, C12sun,N14sun,O16sun};
  int i;

  readSunData();
  fvStat=vfv;
  
  mu=Mcdm*mp_gev/(Mcdm+mp_gev); 
  csSDp=12/M_PI*(pA5*mu)*(pA5*mu)*3.8937966E8*1E-36;  // cm^2
  
  for(sumI=0,i=0;i<17;i++) 
  { double si,FF,fr, cI;
    double vmaxC,vmaxQ;
    MA=A[i]*mp_gev;
    mu=Mcdm*MA/(Mcdm+MA);
    si= (Z[i]*pA0+(A[i]-Z[i])*nA0)*mu;
    si=4/M_PI*si*si*3.8937966E8*1E-36;  // cm^2
    FF=1;
    if(i==0) si+=csSDp;
    if(i<4) {aFraction=sunFractions[i];fr=1;}  else { aFraction=O16sun; fr=A[i]/((double)A[4])*pow(10., P10[i]-P10[4]);}
    muX=(Mcdm-MA)*(Mcdm-MA)/(4*Mcdm*MA);
    vmaxC=(Vesc+Vrot)/(1.E-3*Vlight);
    vmaxQ=vmaxC*vmaxC;
    if(vmaxQ*muX> 2*phiSun[0] ) vmaxQ=2*phiTab[0]/muX;

    if(i==0) FFalpha=0; else FFalpha=Mcdm*MA*pow((0.91*pow(MA,0.3333333) +0.3)*GeVfm,2)/3;
    cI=fr*FF/A[i]*si*simpson(vIntegrand,0,1.E-3*Vlight*sqrt(vmaxQ),1.E-3);     
    sumI+=cI;
//printf("C%d = %E  \n",A[i],rhoDM/Mcdm*cI);
  }
//printf("Sun capture =%E\n",rhoDM*sumI/Mcdm);  
  return  rhoDM*sumI/Mcdm;
}




static double sigmaSun(double pA0, double nA0, double pA5, double nA5, double R3)
{ 
  double si,fr, cI,sumI,csSDp,mu;
  /*                 H    He      C     N     O    Ne     Mg    Si    S     Fe    Na    Al   Cl    Ar     Ca   Cr   Ni */
  int      A[17] ={  1,   4   , 12   ,14   ,   16,   20,   24,   28,   32,   56,   23,   27,   35,   40,   40,  52,  59  };
  int      Z[17] ={  1,   2   ,  6   , 7   ,    8,   10,   12,   14,   16,   26,   11,   13,   17,   18,   20,  24,  28  };
  double P10[17] ={ 12,  10.93,  8.39, 7.78, 8.66, 7.84, 7.53, 7.51, 7.14, 7.45, 6.17, 6.37, 5.50, 6.18, 6.31, 5.64, 6.23};  
  double * sunFractions[5]={Hsun,He4sun, C12sun,N14sun,O16sun};
  int i;
  
  readSunData();
     
  mu=Mcdm*mp_gev/(Mcdm+mp_gev); 
  csSDp=12/M_PI*(pA5*mu)*(pA5*mu)*3.8937966E8*1E-36;  // cm^2

  for(sumI=0,i=0;i<5;i++) 
  {
    MA=A[i]*mp_gev;
    mu=Mcdm*MA/(Mcdm+MA);
    si= (Z[i]*pA0+(A[i]-Z[i])*nA0)*mu;
    si=4/M_PI*si*si*3.8937966E8*1E-36;  // cm^2
    if(i==0) si+=csSDp;
    if(i<4) {aFraction=sunFractions[i];fr=1;}  else { aFraction=O16sun; fr=A[i]/((double)A[4])*pow(10., P10[i]-P10[4]);}

    cI=si/A[i]*(4./3*M_PI)*simpson(sIntegrand,0,R3,1.E-3);   
    sumI+=cI;
  }
  
  cI/=si/A[i]*pow(10,P10[i]);
  
  for(i=5;i<17;i++) 
  { MA=A[i]*mp_gev;
    mu=Mcdm*MA/(Mcdm+MA);
    si= (Z[i]*pA0+(A[i]-Z[i])*nA0)*mu;
    si=4/M_PI*si*si*3.8937966E8*1E-36;  // cm^2
    sumI+=cI*si/A[i]*pow(10,P10[i]);
  } 
  return  sumI;
}

#define RE 6378E5   /* cm */

static double rhoFun(double r) { return polint2(r,nTab,rTab,rhoTab);}

static int readEarthData(void)
{
  static int rdOK=0;
  int i;
  if(!rdOK) 
  { FILE*f;
    char fname[300];
    sprintf(fname,"%s/sources/data_nu/EarthModel.dat",micrO);
    f=fopen(fname,"r"); 
    fscanf(f,"%*[^\n]"); 
    for(i=0;i<SUNPOINTS;i++)
    {
      fscanf(f,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
      rEarth+i, tEarth+i, rhoEarth+i, OEarth+i,  NaEarth+i,  MgEarth+i,   AlEarth+i,  SiEarth+i,  PEarth+i, 
      SEarth+i,   CaEarth+i,  CrEarth+i, FeEarth+i, NiEarth+i);
//    printf("%E %E %E %E %E %E %E %E %E\n", rsun[i], tsun[i], rhosun[i], Hsun[i],He4sun[i], He3sun[i], C12sun[i],N14sun[i],O16sun[i]);
//      printf("frac= %E\n",(1-( Hsun[i]+He4sun[i]+ He3sun[i]+ C12sun[i]+N14sun[i]+O16sun[i]))/O16sun[i] ); 
    }
    fclose(f);
  }   
  rTab= rEarth;
  tTab=tEarth;
  rhoTab=rhoEarth;
  phiTab=phiEarth;
  nTab=EARTHPOINTS;
  Rcm=RE; 
  if(!rdOK) {fillGravPotential(); rdOK=1;} 
  return 0;
}

//double fracFun(double r){ return polint2(r,nTab, rTab,aFraction);}


static double captureEarth(double(*vfv)(double),  double pA0, double nA0, double pA5, double nA5)
{ 
/*            O  Na  Mg   Al  Si  P   S   Ca  Cr, Fe  Ni  Cu  */
  int A[12] ={16, 23 ,24, 27, 28, 30, 32, 40, 52, 56, 58, 64};
  int Z[12] ={ 8, 11 ,12, 13, 14, 15, 16, 20, 24, 26, 28, 29};

  double  sumI,mu;
  double * earthFractions[11]={OEarth,NaEarth,MgEarth,AlEarth,SiEarth,PEarth,SEarth,CaEarth,CrEarth,FeEarth,NiEarth };
  int i;
  
  readEarthData();
  fvStat=vfv;
  
  for(sumI=0,i=0;i<11;i++) 
  { double si,cI;
    double vmaxC,vmaxQ;
    MA=A[i]*mp_gev;
    aFraction=earthFractions[i]; 
    mu=Mcdm*MA/(Mcdm+MA);
    si= (Z[i]*pA0+(A[i]-Z[i])*nA0)*mu;
    si=4/M_PI*si*si*3.8937966E8*1E-36;  // cm^2
    muX=(Mcdm-MA)*(Mcdm-MA)/(4*Mcdm*MA);
    vmaxC=(Vesc+Vrot)/(1.E-3*Vlight);
    vmaxQ=vmaxC*vmaxC;
    if(vmaxQ*muX> 2*phiTab[0] ) vmaxQ=2*phiTab[0]/muX;
    FFalpha=Mcdm*MA*pow((0.91*pow(MA,0.3333333) +0.3)*GeVfm,2)/3; 
    cI=si/A[i]*simpson(vIntegrand,0,1.E-3*Vlight*sqrt(vmaxQ),1.E-3);     
    sumI+=cI;
//printf("C%d = %E  \n",A[i],rhoDM/Mcdm*cI);
  }
  return  rhoDM*sumI/Mcdm;
}



static double sigmaEarth(double pA0, double nA0, double pA5, double nA5, double R3)
{ 
  double si,cI,sumI,mu;
/*               O  Na  Mg   Al  Si  P   S   Ca  Cr, Fe  Ni  Cu  */
  double A[12] ={16, 23 ,24, 27, 28, 30, 32, 40, 52, 56, 58, 64};
  double Z[12] ={ 8, 11 ,12, 13, 14, 15, 16, 20, 24, 26, 28, 29};

  double*earthFractions[11]={OEarth,NaEarth,MgEarth,AlEarth,SiEarth,PEarth,SEarth,CaEarth,CrEarth,FeEarth,NiEarth };
  int i;
  
  readEarthData();
   
  for(sumI=0,i=0;i<11;i++) 
  {
    MA=A[i]*mp_gev;
    mu=Mcdm*MA/(Mcdm+MA);
    si= (Z[i]*pA0+(A[i]-Z[i])*nA0)*mu;
    si=4/M_PI*si*si*3.8937966E8*1E-36;  // cm^2
    aFraction=earthFractions[i];
    cI=si/A[i]*(4./3*M_PI)*simpson(sIntegrand,0,R3,1.E-3);   
    sumI+=cI;
  }  
  return  sumI;
}

double captureAux(int forSun,double(*vfv)(double), double csIp,double csIn,double csDp, double csDn)
{ 
  double pA0, nA0, pA5, nA5;
  double MN=0.939;
  double Mr=MN*Mcdm/(MN+Mcdm);
  double sCoeff= 2*sqrt(3.8937966E8/M_PI)*Mr;
           
  pA0=sqrt(fabs(csIp))/sCoeff;   if(csIp<0) pA0*=-1;
  nA0=sqrt(fabs(csIn))/sCoeff;   if(csIn<0) nA0*=-1;
  pA5=sqrt(fabs(csDp/3))/sCoeff; if(csDp<0) pA5*=-1;
  nA5=sqrt(fabs(csDn/3))/sCoeff; if(csDn<0) nA5*=-1;
  
  if(forSun) return captureSun(vfv,pA0,nA0,pA5,nA5);
       else  return captureEarth(vfv,pA0,nA0,pA5,nA5);               
}

static double  vcs22(numout * cc,int nsub,int * err)
{
   int i;
   double pcm,r;
   double pmass[4], pvect[16];
   double  GG=sqrt(4*M_PI*parton_alpha(3*Mcdm));
   for(i=1;i<=cc->interface->nvar;i++) if(cc->link[i]) 
   cc->interface->va[i]=*(cc->link[i]);

   if( cc->interface->calcFunc()>0 ) {*err=4; return 0;}
   *(cc->interface->gtwidth)=0;
   *(cc->interface->twidth)=0;
   *(cc->interface->gswidth)=0;
   for(i=0;i<4;i++) cc->interface->pinf(nsub,1+i,pmass+i,NULL);
   *err=0;
   if(pmass[0]+pmass[1] <= pmass[2]+pmass[3]) return 0;
   for(i=0;i<16;i++) pvect[i]=0;

   pcm= decayPcm(pmass[0]+pmass[1],pmass[2],pmass[3]);
   for(i=0;i<2; i++) pvect[4*i]=pmass[i];
   for(i=2;i<4; i++) pvect[4*i]=sqrt(pmass[i]*pmass[i] +pcm*pcm);
   pvect[8+3]=pcm;
   pvect[12+3]=-pcm;
   r=cc->interface->sqme(nsub,GG,pvect,err);
   return 3.8937966E8*r*pcm/(16*M_PI*pmass[0]*pmass[1]*(pmass[0]+pmass[1]));
}



/* ================ Basic Spectra  =====   hep-ph/0506298v5 ============== */


static double massTabSun[12]={10, 30, 50, 70, 90, 100, 200, 300, 500, 700, 900, 1000};
static double xTab[100]={0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.2,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,0.3,0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.39,0.4,0.41,0.42,0.43,0.44,0.45,0.46,0.47,0.48,0.49,0.5,0.51,0.52,0.53,0.54,0.55,0.56,0.57,0.58,0.59,0.6,0.61,0.62,0.63,0.64,0.65,0.66,0.67,0.68,0.69,0.7,0.71,0.72,0.73,0.74,0.75,0.76,0.77,0.78,0.79,0.8,0.81,0.82,0.83,0.84,0.85,0.86,0.87,0.88,0.89,0.9,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99,1.};
static double nuTabSun[12][2][8][100];
static char * chan[8]={"nu","b","tau","c","q","t","Z","W"};
static double nuTabEarth[14][2][8][100];
static double massTabEarth[14]= {10, 30, 50, 70, 90, 100, 150, 200, 250, 300, 500, 700, 900, 1000};

static int rdSunOk=0,rdEarthOk=0;

static int   nuSunRead(void)
{ int i,j,k,l,nLine=2;
  double x,mass;
  char *fname;
  const char* fnameL[2]={"Sun_numu_EVOL.dat","Sun_numub_EVOL.dat"};
  
  fname=malloc(strlen(micrO)+100);
  for(l=0;l<2;l++)
  {
     FILE*f;
     sprintf(fname,"%s/sources/data_nu/%s",micrO,fnameL[l]);
     f=fopen(fname,"r"); 
     if(!f) return -(l+1);
     fscanf(f,"%*[^\n]");   
    
     for(i=0;i<12;i++) for(j=0;j<100;j++) 
     { fscanf(f," %lf %lf ",&mass, &x);
       if(mass!=massTabSun[i]) { printf(" wrong mass %E !=%E  line %d\n",mass,massTabSun[i],nLine); return 1;}
       if(x!=xTab[j])          { printf(" wrong x  %E !=%E line=%d  \n",x,xTab[j],nLine); return 2;}
       for(k=0;k<8;k++) fscanf(f," %lf ",nuTabSun[i][l][k]+j);
       nLine++;            
     }
     fclose(f);
  }
  free(fname);
  rdSunOk=1;
  return 0;
}  

static int   nuEarthRead(void)
{ int i,j,k,l,nLine=2;
  double x,mass;
  char *fname=malloc(strlen(micrO)+100);
  const char * fnameL[2]={"Earth_numu_EVOL.dat","Earth_numub_EVOL.dat"};
  for(l=0;l<2;l++)
  {
     FILE*f;
     sprintf(fname,"%s/sources/data_nu/%s",micrO,fnameL[l]); 
     f=fopen(fname,"r");
     if(!f) return -(l+1);
     fscanf(f,"%*[^\n]");   
    
     for(i=0;i<14;i++) for(j=0;j<100;j++) 
     { fscanf(f," %lf %lf ",&mass, &x);
       if(mass!=massTabEarth[i]) { printf(" wrong mass %E !=%E  line %d\n",mass,massTabEarth[i],nLine); return 1;}
       if(x!=xTab[j])       { printf(" wrong x  %E !=%E line=%d  \n",x,xTab[j],nLine); return 2;}
       nuTabEarth[i][l][0][j]=0;
       for(k=1;k<8;k++) fscanf(f," %lf ",nuTabEarth[i][l][k]+j);
       nLine++;            
     }
     fclose(f);
  }
  free(fname);
  rdEarthOk=1;
  return 0;  
}


static void mInterpSun(double Nmass,  int  CHin, int  CHout, double*tab)
{  
   int l,i0;
   double c0,c1;
   double *p0,*p1;
   if(!rdSunOk) nuSunRead();
   for(i0=0; i0<12 && Nmass>=massTabSun[i0] ;i0++);
   if(i0) i0--;
   p0=nuTabSun[i0  ][CHout][CHin];
   p1=nuTabSun[i0+1][CHout][CHin];
   if(i0==12-1) for(l=0;l<100;l++) tab[l]= p0[l];
   else
   {   
     c1=(Nmass*Nmass -massTabSun[i0]*massTabSun[i0])/(massTabSun[i0+1]*massTabSun[i0+1] - massTabSun[i0]*massTabSun[i0]);
     c0=1-c1;
     for(l=0;l<100;l++) tab[l]= c0*p0[l]+c1*p1[l];
   }
}   


static void mInterpEarth(double Nmass,  int  CHin,int  CHout, double*tab)
{  
   int l,i0;
   double c0,c1;
   double *p0,*p1;
   if(!rdEarthOk) nuEarthRead();
   for(i0=0; i0<14 && Nmass>=massTabEarth[i0] ;i0++);
   if(i0) i0--;
   p0=nuTabEarth[i0  ][CHout][CHin];
   p1=nuTabEarth[i0+1][CHout][CHin];
   if(i0==14-1) for(l=0;l<100;l++) tab[l]= p0[l];
   else
   {   
     c1=(Nmass*Nmass -massTabEarth[i0]*massTabEarth[i0])/(massTabEarth[i0+1]*massTabEarth[i0+1] - massTabEarth[i0]*massTabEarth[i0]);
     c0=1-c1;
     for(l=0;l<100;l++) tab[l]= c0*p0[l]+c1*p1[l];
   }
}   
  
                                                     
static int basicNuSpectra_(int forSun, double M, int pdgN, int outN, double * tab)
{ 
  int inP,i,j;
  int N=abs(pdgN);
  double tab100[100];
  double c=1;
  
  for(i=0;i<NZ;i++) tab[i]=0;
  
//  if(N==12 ||N==14 ||N==16) {  if( pdgN*outN<0) return 0; else c=1;}
  
  switch(N)
  {
    case 1: case 2: 
    case 3: case 21: inP=4; break;
    case 4:          inP=3; break;
    case 5:          inP=1; break; 
    case 6:          inP=5; break;
 
    case 15:         inP=2; break;  /*l  */ 
    case 23:         inP=6; break;  /*z  */ 
    case 24:         inP=7; break;  /*w  */
    case 12:
    case 14:
    case 16:  if(forSun) inP=0; else 
              {  if((pdgN==14 && outN>0)||(pdgN==-14 && outN<0)) tab[0]=2/(Zi(0)-Zi(1)); 
                 return 0; 
              }     /*nu */
              break;  
    default:                       return 1;
  }  
  
  if(forSun) mInterpSun(M/2,inP,-(outN-1)/2,tab100); 
    else     mInterpEarth(M/2,inP,-(outN-1)/2,tab100);

  for(i=0;i<NZ;i++)
  { double x=exp(Zi(i));
    tab[i]+=c*x*polint3(x,100, xTab ,tab100);  
  }

  return 0;
}

int basicNuSpectra(int forSun, int pdgN, int outN, double * tab)
{ return  basicNuSpectra_(forSun,2*Mcdm, pdgN, outN, tab);}


/*  ===================  Spectra ========================== */


static void getSpectrum(int forSun, double M, double m1,double m2,char*n1,char*n2, int N1, int N2,int outP, double *tab)
{
  int i,k;
  char* nn[2];
  int pdg[2];
  double mm[2],E[2],p2;

  
  for(i=0;i<NZ;i++) tab[i]=0; 

  if(abs(N1)==abs(N2)) switch(abs(N1))
  { case 22: case 11: case 13: return;}
     
  if(N1+N2==0  || (N1==21 && N2==21) || (N1==23 && N2==23) )  { if(basicNuSpectra_(forSun, M, N1, outP,tab)==0) return;} 

  nn[0]=n1;  nn[1]=n2;
  mm[0]=m1;  mm[1]=m2;  
  pdg[0]=N1; pdg[1]=N2;          
  
  
  if(M>m1+m2) p2=sqrt((M*M-(m1+m2)*(m1+m2))*(M*M-(m1-m2)*(m1-m2)))/(2*M);
  else 
  { p2=0; 
    if(abs(N1)==abs(N2)) { mm[0]=M/2; mm[1]=M/2;} else
    {
         if(N1==23 || abs(N1)==24)   mm[0]=M-m2; else mm[1]=M-m1;
    }   
  }
    
  E[0]=sqrt(mm[0]*mm[0]+p2*p2);
  E[1]=sqrt(mm[1]*mm[1]+p2*p2);
    
  for(k=0;k<2;k++)
  {   double dY;
      double tabAux[NZ];
      
      if(abs(pdg[k])==11 || abs(pdg[k])==13 || abs(pdg[k])==22) continue;
      if(basicNuSpectra_(forSun,2*E[k],pdg[k], outP,tabAux)==0)
      { double kf; 
        dY=log(M/E[k]/2);
        switch(abs(pdg[k]))
        { case 12: case 14: case 16: 
          if(pdg[k]*outP>1) kf=1; else kf=0; break;
          default: kf=0.5;
        } 
        boost(dY, M/2, E[k], 0., tabAux);
        for(i=0;i<NZ;i++)tab[i]+=kf*tabAux[i]; 
      }
      else
      { double w=0;
        numout * d2Proc;
        int l; 
        char* n[4];
        double m[4];
        double tab_p[NZ];
        char process[40],plib[40];
        int ntot;

        if(mm[k]==0) { fprintf(stderr,"Can not hadronize BSM zero mass %s\n",nn[k]); continue;}
                                         
        strcpy(plib,"2width_");
        sprintf(process,"%s->2*x",nn[k]);
        pname2lib(nn[k],plib+7);
        for(i=0;i<NZ;i++) tabAux[i]=0;
                     
        d2Proc=getMEcode(0,ForceUG,process,NULL,NULL,plib);
        if(!d2Proc) { fprintf(stderr,"Can not find decay modes for  mass %s\n",nn[k]); continue; }
        procInfo1(d2Proc,&ntot,NULL,NULL);
        { double  Qstat;
          if(Qaddress) { Qstat=*Qaddress; setQforParticle(Qaddress,nn[k]); }  
          for(l=1;l<=ntot ;l++)
          {    
            double wP=pWidth2(d2Proc,l);
            
            if(wP>0)
            { int N2=d2Proc->interface->pinfAux(l,2,NULL,NULL,NULL);
              int N3=d2Proc->interface->pinfAux(l,3,NULL,NULL,NULL); 
              procInfo2(d2Proc,l,n,m);    
              getSpectrum(forSun,m[0],m[1],m[2],n[1],n[2],N2,N3,outP, tab_p);
              for(i=0;i<NZ;i++) tabAux[i]+=wP*tab_p[i];
              w+=wP;
            }
          }
          if(Qaddress){ *Qaddress=Qstat; calcMainFunc();}
        }
         
        if(w==0) { if(abs(pdg[k])!= abs(pNum(CDM)))   fprintf(stderr,"Can't find decays for  %s\n",nn[k]);
                   continue;
                 }
        dY=acosh(E[k]/mm[k]);
        boost(dY, M/2, mm[k], 0., tabAux);
        for(i=0;i<NZ;i++)tab[i]+=tabAux[i]/w;
      }
    } 
}



static double calcSpectrum0(char *name1,char*name2, int forSun,   double *Spectranu, double *SpectraNu, double *alpha)
{
  int i,k;
  double vcsSum=0,vcsSum1=0; 
  int ntot,err;
  double * v_cs;
  
  char name1L[10],name2L[10], lib[20],process[400];
  numout * libPtr;
  
  for(i=0;i<NZ;i++) Spectranu[i]=SpectraNu[i]=0;  

  pname2lib(name1,name1L);
  pname2lib(name2,name2L);
  sprintf(lib,"omg_%s%s",name1L,name2L);
  sprintf(process,"%s,%s->AllEven,1*x{%s",name1,name2,EvenParticles());
// Warning!!   in should be done in the same manner as annihilation libraries for Omega
  libPtr=getMEcode(0,ForceUG,process,NULL,NULL,lib);
  
  
  if(!libPtr) return 0;
  passParameters(libPtr);
  procInfo1(libPtr,&ntot,NULL,NULL); 
  
  v_cs=malloc(sizeof(double)*ntot);
  (*libPtr->interface->twidth)=0;
  
  for(k=0;k<ntot;k++)
  { double m[4];
    char *N[4];
    int pdg[4];
    int l,l_;
    double br,wV;
    
    for(i=0;i<4;i++) N[i]=libPtr->interface->pinf(k+1,i+1,m+i,pdg+i);
    cc23=NULL;
    v_cs[k]=0;

    if(VZdecay||VWdecay)
    {  int nVV;
       int vd[4]={0,0,0,0};
       for(l=2;l<4;l++) if((pdg[l]==23&&VZdecay) || (abs(pdg[l])==24&&VWdecay)) vd[l]=1;
            
       for(l=2;l<4;l++) if(vd[l]) break;
       if(l<3)
       {  l_=5-l; 
          if(vd[l_])
          { nVV=2;
            if(m[l_]>m[l]) { l=l_; l_=5-l;}
          } else nVV=1;
          
          if(m[0]+m[1] >  m[l_] +20  && m[0] + m[1] <  m[2]+m[3] + 4*nVV)
           cc23=xVtoxll(2,2,N,pdg,l,&wV,&br);                
       }
    }
    if(cc23)
    { int i3W;  
      double  r,m1,v0=0.001;
      for(i3W=2;i3W<5;i3W++) if(strcmp(cc23->interface->pinf(1,i3W+1,NULL,NULL),N[l_])==0) break;
      r=v0*cs23(cc23,1,v0*Mcdm/2,i3W)/br;
       
      if(pdg[l_]==23 || abs(pdg[l_])==24)
      { double wV2;
        
        wV2=pWidth(N[l_],NULL);
        r*=decayPcmW(2*Mcdm,m[l],m[l_],wV,wV2,0)/decayPcmW(2*Mcdm,m[l],m[l_],wV,0,0);
        if(pdg[l]==pdg[l_]) r/=2;
      }
      v_cs[k]=r;
      vcsSum+=r;                     
    }
    else if((m[2]+m[3])<m[0]+m[1])
    { 
#ifdef V0    
      v_cs[k]=V0*cs22(libPtr,k+1,V0*m[0]/2,-1.,1.,&err);
#else 
      v_cs[k]= vcs22(libPtr,k+1,&err);
#endif 
      if(v_cs[k]<0) v_cs[k]=0; 
      vcsSum+=v_cs[k];
    } else v_cs[k]=-1;
  }
   
  for(k=0;k<ntot ;k++) if(v_cs[k]>0)
  { char * N[4];
    double m[4];
    int l,pdg[2];
    int PlusAok=0;
    double tab2[NZ];

    procInfo2(libPtr,k+1,N,m);
    for(l=0;l<2;l++)  pdg[l]=qNumbers(N[2+l],NULL,NULL,NULL);
    if(N[2][0]=='~' || N[3][0]=='~') vcsSum1+=v_cs[k]; 
#ifdef PRINT
       { char txt[100];
         sprintf(txt,"%s,%s -> %s %s", N[0],N[1],N[2],N[3]);
         printf("  %-20.20s  %.2E\n",txt,v_cs[k]*2.9979E-26);
       }
#endif       
    getSpectrum(forSun, pMass(name1)+pMass(name2), m[2], m[3],N[2],N[3],pdg[0],pdg[1], 1,tab2);
    for(i=0;i<NZ;i++) Spectranu[i]+=tab2[i]*v_cs[k]/vcsSum;
    getSpectrum(forSun, pMass(name1)+pMass(name2), m[2], m[3],N[2],N[3],pdg[0],pdg[1],-1,tab2);
    for(i=0;i<NZ;i++) SpectraNu[i]+=tab2[i]*v_cs[k]/vcsSum;  
  } 
  free(v_cs);
  if(alpha) { if(vcsSum>0) *alpha=vcsSum1/vcsSum; else *alpha=0;}
  return  vcsSum*2.9979E-26;
}

/* ========  Sun[Earth] neutrino fluxes ============= */

#define Gconst (0.7426E-30) /* m/g */
#define KelvinEv (8.61734E-05) /* ev */
#define Etime  (1.5E17)     /* time of existence of Sun and Earth in seconds */ 

static double C[2],An[2],Ev[2],Alpha;

static void deriv1(double t,double*n,double*dn) { dn[0]=C[0]-An[0]*n[0]*n[0]-Ev[0]*n[0];

//printf(" C=%E An=%E Ev=%E\n", C[0], An[0]*n[0]*n[0], Ev[0]*n[0]);
} 
  
static void deriv2(double t,double*n,double*dn) 
{ dn[0]=C[0]+ 0.5*Alpha*An[1]*n[1]*n[1] - An[1]*n[0]*n[0] - An[0]*n[0]*n[1] -Ev[0]*n[0]; 
  dn[1]=C[1]+ 0.5*Alpha*An[1]*n[0]*n[0] - An[1]*n[1]*n[1] - An[0]*n[0]*n[1] -Ev[1]*n[1];
}

int neutrinoFlux(double (* fvf)(double), int forSun, double* nu, double * Nu)
{
  int i,n,err;
  double vcs0,vcs1;

  double nu_[NZ],Nu_[NZ];
  char *name, *aname;
  double pA0[2],pA5[2],nA0[2],nA5[2];
  double R,Prop,Cr0,Cr1,Dv; 
  double Veff;
  double rho,T;

  for(i=0;i<NZ;i++) nu[i]=Nu[i]=0;

  err=nucleonAmplitudes(FeScLoop,pA0,pA5,nA0,nA5);
  if(err) return err;
  Cr0=forSun? captureSun(fvf,pA0[0],nA0[0],pA5[0],nA5[0]):captureEarth(fvf, pA0[0],nA0[0],pA5[0],nA5[0]); 
  
  if(pA0[0]==pA0[1] && nA0[0]==nA0[1] &&pA5[0]==pA5[1]&&nA5[0]==nA5[1]) Cr1=Cr0; else
  Cr1=forSun? captureSun(fvf,pA0[1],nA0[1],pA5[1],nA5[1]):captureEarth(fvf,pA0[1],nA0[1],pA5[1],nA5[1]);  
                         
  name=CDM;
  aname=pdg2name(-pNum(CDM));
  if(!aname) aname=name;
  if(forSun) R=150E6;  else R=6378.1; /* Distance to Sun/Earth in [km] */  
  Prop=31556925.2/(4*M_PI*R*R);       /* for Year*km^2 */
    
  vcs0= calcSpectrum0(name,aname,forSun, nu,Nu,NULL);
  { 
     double r_,v_,r095,S,ph,Veff1;
     for(i=0,r_=0;i<10;i++) 
     { T=polint2(r_/Rcm,nTab,rTab,tTab)*KelvinEv*1E-9;
       rho= polint2(r_/Rcm,nTab,rTab,rhoTab);
       r_= sqrt(6*T/(Gconst*100*rho*Mcdm))/M_PI;
       if(r_>Rcm) r_=Rcm;
     }

     rho=polint2(r_/Rcm,nTab,rTab,rhoTab);
     T=polint2(r_/Rcm,nTab,rTab,tTab);
     r095=polint2(T*0.95,nTab,tTab,rTab);
     if(r095>1) r095=1;
     T*=KelvinEv*1E-9;
     r095*=Rcm;
     Veff1=pow(r_*M_PI,3)/8;
     if(forSun) S=sigmaSun(pA0[0],nA0[0],pA5[0],nA5[0],r095*r095*r095);
     else       S=sigmaEarth(pA0[0],nA0[0],pA5[0],nA5[0],r095*r095*r095);
     v_=sqrt(8*T/(M_PI*Mcdm));
     ph=phiTab[0]*Mcdm/T;
     ph=ph*exp(-ph);
     
     Ev[0]=v_*S/Veff1*ph*Vlight*100;
     if(strcmp(name,aname))
     { if(forSun)S=sigmaSun(pA0[1],nA0[1],pA5[1],nA5[1],r095);
       else S=sigmaEarth(pA0[1],nA0[1],pA5[1],nA5[1],r095);
       Ev[1]=v_*S/Veff1*ph*Vlight*100; 
     } 
     Veff=pow(Gconst*100*rho*Mcdm/T/3,-1.5);
  } 

  if(strcmp(name,aname))
  { double G01,G00,G11;
    double N[2]={0,0};
    int err;
    vcs1=calcSpectrum0(name,name,forSun, nu_,Nu_,&Alpha);
    
    C[0]=Cr0*(1+dmAsymm)/2;                             
    C[1]=Cr1*(1-dmAsymm)/2;       
       
    An[0]=vcs0/Veff;
    An[1]=vcs1/Veff;
    err=odeint(N,2, 0 ,Etime , 1.E-3,  Etime/10, deriv2);
//printf("err=%d N={%E,%E}\n",err, N[0],N[1]);
    G00=0.5*An[1]*N[0]*N[0];
    G11=0.5*An[1]*N[1]*N[1];
    G01=    An[0]*N[0]*N[1];    
    { /* symbolic solution for large Etime */
      double alpha,beta,x,G01_,G00_,G11_;        
      alpha=vcs1/vcs0;
      beta= Cr0/Cr1;
      x=(beta-1 + sqrt((beta-1)*(beta-1) + 4*beta*alpha*alpha))/(2*alpha);
      /* x= rho_particle/rho_antiparticle  */ 
      G01_=Cr0/(1+alpha*x);
      G00_=0.5*G01_*alpha*x;
      G11_=0.5*G01_*alpha/x;
//      printf("x=%E G00 = %E/%E  G01 = %E/%E G11 = %E/%E\n",x, G00,G00_,G01,G01_,G11,G11_);
    }
    for(i=0;i<NZ;i++) 
    { nu[i]=Prop*(G01*nu[i]+G00*nu_[i]+G11*Nu_[i]); 
      Nu[i]=Prop*(G01*Nu[i]+G00*Nu_[i]+G11*nu_[i]);
    }  
    vcs1=calcSpectrum0(aname,aname,forSun, nu_,Nu_,NULL);     
    for(i=0;i<NZ;i++)
    { nu[i]+=Prop*(G11*nu_[i]);
      Nu[i]+=Prop*(G11*Nu_[i]);
    }              

  } else   
  {   
    int err;  
    double N=0;
    double rf;
    C[0]=Cr0;
    An[0]=vcs0/Veff;
    err=odeint(&N,1, 0 ,Etime , 1.E-3,  Etime/10, deriv1);
    rf=0.5*An[0]*N*N*Prop;

//printf("Rate Factor = %E(num.sol), =%E(formula)\n",rf,   0.5*Cr0*pow(tanh(Etime*sqrt(Cr0*vcs0/Veff)),2)*Prop);
    for(i=0;i<NZ;i++) { nu[i]*=rf;  Nu[i]*=rf;}
  }
  return 0;
}

/* ==================  Muon flux =========  0906.4364  =============  */

static double *SpN_stat=NULL;
static double tabNuSpectrum(double E) {return SpectdNdE(E,SpN_stat);}
static double (*nuSpectrum)(double)=tabNuSpectrum;
static double Enu_stat,Emu_stat,alpha_stat,beta_stat;


static int inPr=1; /* proton;   -1 for neutron */
static int inNu=1; /* neutrino; -1 for anti-neutrino */


static double dSigmadE_nu2mu(double Enu,double Emu)   /*result in   1/cm^2/GeV */ 
{  double GF=1.16637E-5;  /* GeV^(−2) */
   double GeVcm=0.50677E14;
   double C=2*GF*GF*mp_gev/M_PI/GeVcm/GeVcm;
   double a=0,b=0;
   switch (inNu)
   { case  1:switch (inPr)
             { 
               case  1: a=0.15,b=0.04; break;
               case -1: a=0.25,b=0.06; break; 
             }
             break;  
     case -1:switch (inPr)
             { 
               case  1: a=0.06,b=0.25; break;
               case -1: a=0.04,b=0.15; break; 
             }
             
             break;
   }
   return C*(a+b*Emu*Emu/(Enu*Enu));       
} 

static double integrandX(double x)
{  double g=alpha_stat/beta_stat;  
   double ex=exp(x*beta_stat);
   double q= (0.10566 /* muon mass*/)/(65865.4 /* (muon_life_time)*c in cm */)/alpha_stat;
   double Emu_prim=(Emu_stat+g)*ex-g;
   double Psurv=pow(Emu_stat*(Emu_prim+g)/(Emu_prim*(Emu_stat+g)),q);
   return  dSigmadE_nu2mu(Enu_stat, Emu_prim )*ex*Psurv; 
} 

static double integrandEnuUpward(double e)
{ double E;
  if(e==0) return 0;
  E=1/e;
  Enu_stat=E;
  double g=alpha_stat/beta_stat;
  double xMax=log((Enu_stat+g)/(Emu_stat+g))/beta_stat;
  return  E*E*nuSpectrum(E)*simpson(integrandX,0,xMax,1.E-4);
}

void muonUpward(double*nu, double*Nu, double*mu)
{
   int i,k,l;
   double C;
   double  NA=6.002141E23 /* mol-1*/;
   double Nrate[2]={0.5,0.5}; /* proton , neutron */
   double*Sp[2]={nu,Nu};
   double rho=2.6;
   
   nuSpectrum=tabNuSpectrum;  
   alpha_stat=0.002*rho;
   beta_stat=3.0E-6*rho;

   mu[0]=0;
   for(i=1;i<NZ;i++) 
   {  Emu_stat=Mcdm*exp(Zi(i));
      mu[i]=0;
      if(Emu_stat>0.01) for(k=0;k<2;k++) for(l=0;l<2;l++)
      {
         inNu=1-2*k;
         if(Sp[k])
         {
           SpN_stat=Sp[k];
           inPr=-1+2*l;
           mu[i]+=simpson(integrandEnuUpward,1/Mcdm,1/Emu_stat,1.E-4)*Nrate[l];          
         }  
      }   
      mu[i]*=rho*NA*Emu_stat;      
   }
}

static double integrandEnuContained(double e) 
{  double E; 
   if(e==0) return 0; 
   E=1/e; return  E*E*nuSpectrum(E)*dSigmadE_nu2mu(E,Emu_stat) ; 
}

void muonContained(double*nu,double*Nu,double rho, double*mu)
{
  int i,k,l;
  double C;
  double  NA=6.002141E23 /* mol-1*/;
  double Nrate[2]={0.5,0.5}; /* proton , neutron */
  double*Sp[2]={nu,Nu};
  
  nuSpectrum=tabNuSpectrum; 
  mu[0]=0; 

  for(i=1;i<NZ;i++)
  {  Emu_stat=Mcdm*exp(Zi(i));
     mu[i]=0;
     if(Emu_stat>0.01) for(k=0;k<2;k++) for(l=0;l<2;l++)
     { inNu=1-2*k;
       if(Sp[k])
       {
         SpN_stat=Sp[k];
         inPr=-1+2*l;
         mu[i]+=simpson(integrandEnuContained,1/Mcdm, 1/Emu_stat,1.E-4)*Nrate[l];
       }  
     }  
     mu[i]*=rho*NA*Emu_stat*1E5;
  }                  
}

// Background 

static double cosFi_stat;

static double  ATMdNudE(double E) // for NuBar *1.35/1.95 
{
  int i;
  const double gamma=1.74, a=0.018, b=0.024,c=0.0069, e=0.00139;
  return 1.95E17*pow(E,-gamma-1)*(a/(1+b*E*cosFi_stat)+c/(1+e*E*cosFi_stat));
} 


double  ATMmuonUpward(double cosFi, double E)
{ 
   int k,l;
   double  NA=6.002141E23 /* mol-1*/;
   double Nrate[2]={0.5,0.5}; /* proton , neutron */
   double rho=2.6;
   double mu=0;
     
   alpha_stat=0.002*rho;
   beta_stat=3.0E-6*rho;
   Emu_stat=E;
   cosFi_stat=cosFi;
   nuSpectrum=ATMdNudE;
   
   for(k=0;k<2;k++) for(l=0;l<2;l++)
   {
      inNu= 1-2*k;
      inPr=-1+2*l;
      mu+= (k==0?1: 1.35/1.95)*simpson(integrandEnuUpward,0,1/E,1.E-4)*Nrate[l];
   }   
   return  mu*rho*NA;
}

double  ATMmuonContained(double cosFi, double E,double rho)
{ 
   int k,l;
   double  NA=6.002141E23 /* mol-1*/;
   double Nrate[2]={0.5,0.5}; /* proton , neutron */
   double mu=0;
   
   alpha_stat=0.002*rho;
   beta_stat=3.0E-6*rho;
   Emu_stat=E;
   cosFi_stat=cosFi;
   nuSpectrum=ATMdNudE;
   
   for(k=0;k<2;k++) for(l=0;l<2;l++)
   {
      inNu= 1-2*k;
      inPr=-1+2*l;
      mu+= (k==0?1: 1.35/1.95)*simpson(integrandEnuContained,0,1/E,1.E-4)*Nrate[l];
   }   
   return  mu*rho*NA*1.E5;      
}
