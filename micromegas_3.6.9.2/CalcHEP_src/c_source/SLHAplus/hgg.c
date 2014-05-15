#include"SLHAplus.h"

//#define HDECAY

static double   fiRe(double tau)
{
  double x;
  if(tau<1)
  { x=asin(sqrt(tau));
    return x*x;
  }else if(tau==1) return 0;
  else
  {
    x=sqrt(1-1/tau);   
    x=log((1+x)/(1-x));
    return -0.25*(x*x - M_PI*M_PI);
  }
}

static double fiIm(double tau)
{
  double x;
  if(tau<=1) return 0;
  else
  {
     x=sqrt(1-1/tau);   
     x=log((1+x)/(1-x));
     return 0.5*x*M_PI;
  }
}

static double HggFr(double tau) { return  2*(tau+(tau-1)*fiRe(tau))/(tau*tau); }
static double HggFi(double tau) { return  2*(tau-1)*fiIm(tau)/(tau*tau); }
static double HggVr(double tau) { return -2*(2*tau*tau + 3*tau + 3*(2*tau-1)*fiRe(tau))/(2*tau*tau); }
static double HggVi(double tau) { return -2*(3*(2*tau-1)*fiIm(tau))/(2*tau*tau); }
static double HggAr(double tau) { return  2*fiRe(tau)/tau;}
static double HggAi(double tau) { return  2*fiIm(tau)/tau;}
static double HggSr(double tau) { return -(tau-fiRe(tau))/(tau*tau);}
static double HggSi(double tau) { return  fiIm(tau)/(tau*tau);}


double complex HggF(double tau) { return HggFr(tau)+I*HggFi(tau);}
double complex HggV(double tau) { return HggVr(tau)+I*HggVi(tau);}
double complex HggA(double tau) { return HggAr(tau)+I*HggAi(tau);}
double complex HggS(double tau) { return HggSr(tau)+I*HggSi(tau);}


#ifdef HDECAY
#include<complex.h>

extern void ckofq_hdec__(double complex * res,  int * icase, double * rho);
extern void ckofsq_hdec__(double complex * res,  double * rho);

double complex  hdHgam1Sc(double x) 
{ double complex res;
  double rho=x*4;
  ckofsq_hdec__(&res,&rho);
  return res;
}

double complex  hdHgam1Fc(double x) 
{ double complex res;
  double rho=x*4;
  int icase=0;
  ckofq_hdec__(&res,&icase,&rho);
  return res;
}

double complex  hdHgam1Ac(double x) 
{ double complex res;
  double rho=x*4;
  int icase=1;
  ckofq_hdec__(&res,&icase,&rho);
  return res;
}

double   hdHgam1Fr(double x){ return creal(hdHgam1Fc(x));}
double   hdHgam1Ar(double x){ return creal(hdHgam1Ac(x));}
double   hdHgam1Sr(double x){ return creal(hdHgam1Sc(x));}
double   hdHgam1Fi(double x){ return cimag(hdHgam1Fc(x));}
double   hdHgam1Ai(double x){ return cimag(hdHgam1Ac(x));}
double   hdHgam1Si(double x){ return cimag(hdHgam1Sc(x));}


double  alphas_hdec__(double *Q, int *N) { return alphaQCD(*Q);}

#endif


#define mNN 71 
#define maxTauTab 9.8656E+01

static double mX[mNN]={  0.0000E+00, 5.0000E-02, 1.0000E-01, 1.5000E-01, 2.0000E-01, 2.5000E-01, 3.0000E-01, 3.5000E-01, 4.0000E-01, 4.5000E-01, 5.0000E-01, 5.5000E-01, 6.0000E-01, 6.5000E-01, 7.0000E-01, 7.5000E-01, 8.0000E-01, 8.5000E-01, 9.0000E-01, 9.5000E-01, 1.0000E+00, 1.0500E+00, 1.1000E+00, 1.1500E+00, 1.2000E+00, 1.2500E+00, 1.3000E+00, 1.3500E+00, 1.4000E+00, 1.4500E+00, 1.5000E+00, 1.5500E+00, 1.6000E+00, 1.6500E+00, 1.7000E+00, 1.7500E+00, 1.8000E+00, 1.8500E+00, 1.9000E+00, 1.9500E+00, 2.0000E+00, 2.0500E+00, 2.1000E+00, 2.1500E+00, 2.2000E+00, 2.2500E+00, 2.3000E+00, 2.3500E+00, 2.4000E+00, 2.4500E+00, 2.5000E+00, 2.5500E+00, 2.6000E+00, 2.6500E+00, 2.7000E+00, 2.7500E+00, 2.8000E+00, 2.8500E+00, 2.9000E+00, 2.9500E+00, 3.0000E+00, 3.0500E+00, 3.1000E+00, 3.1500E+00, 3.2000E+00, 3.2500E+00, 3.3000E+00, 3.3500E+00, 3.4000E+00, 3.4500E+00, 3.5000E+00};
//            tau[mNN]={ 0.0000E+00, 2.2622E-01, 4.0951E-01, 5.5629E-01, 6.7232E-01, 7.6270E-01, 8.3193E-01, 8.8397E-01, 9.2224E-01, 9.4967E-01, 9.6875E-01, 9.8155E-01, 9.8976E-01, 9.9475E-01, 9.9757E-01, 9.9902E-01, 9.9968E-01, 9.9992E-01, 9.9999E-01, 1.0000E+00, 1.0000E+00, 1.0000E+00, 1.0000E+00, 1.0001E+00, 1.0003E+00, 1.0010E+00, 1.0024E+00, 1.0053E+00, 1.0102E+00, 1.0185E+00, 1.0312E+00, 1.0503E+00, 1.0778E+00, 1.1160E+00, 1.1681E+00, 1.2373E+00, 1.3277E+00, 1.4437E+00, 1.5905E+00, 1.7738E+00, 2.0000E+00, 2.2763E+00, 2.6105E+00, 3.0114E+00, 3.4883E+00, 4.0518E+00, 4.7129E+00, 5.4840E+00, 6.3782E+00, 7.4097E+00, 8.5938E+00, 9.9466E+00, 1.1486E+01, 1.3230E+01, 1.5199E+01, 1.7413E+01, 1.9896E+01, 2.2670E+01, 2.5761E+01, 2.9195E+01, 3.3000E+01, 3.7205E+01, 4.1841E+01, 4.6940E+01, 5.2536E+01, 5.8665E+01, 6.5363E+01, 7.2670E+01, 8.0626E+01, 8.9274E+01, 9.8656E+01};
static double mFr[mNN]={-1.0000E+00,-5.3444E-01,-2.0534E-01, 8.0299E-02, 3.3935E-01, 5.7820E-01, 7.9923E-01, 1.0028E+00, 1.1881E+00, 1.3536E+00, 1.4978E+00, 1.6191E+00, 1.7169E+00, 1.7914E+00, 1.8442E+00, 1.8782E+00, 1.8975E+00, 1.9067E+00, 1.9087E+00, 1.9078E+00, 1.9075E+00, 1.9088E+00, 1.9131E+00, 1.9170E+00, 1.9233E+00, 1.9419E+00, 1.9733E+00, 2.0195E+00, 2.0781E+00, 2.1401E+00, 2.1849E+00, 2.2174E+00, 2.1993E+00, 2.1260E+00, 1.9954E+00, 1.8128E+00, 1.5924E+00, 1.3474E+00, 1.0954E+00, 8.4661E-01, 6.0795E-01, 3.8943E-01, 1.9104E-01, 1.5805E-02,-1.3671E-01,-2.6851E-01,-3.8171E-01,-4.7869E-01,-5.5836E-01,-6.2561E-01,-6.8028E-01,-7.2488E-01,-7.6087E-01,-7.8994E-01,-8.1277E-01,-8.2981E-01,-8.4217E-01,-8.4947E-01,-8.5388E-01,-8.5845E-01,-8.5956E-01,-8.5203E-01,-8.5102E-01,-8.4926E-01,-8.4304E-01,-8.3484E-01,-8.2734E-01,-8.1977E-01,-8.1058E-01,-8.0218E-01,-7.8552E-01};
static double mAr[mNN]={ 0.0000E+00, 7.2793E-01, 1.3009E+00, 1.8526E+00, 2.4126E+00, 2.9981E+00, 3.6229E+00, 4.3002E+00, 5.0449E+00, 5.8742E+00, 6.8092E+00, 7.8776E+00, 9.1154E+00, 1.0572E+01, 1.2319E+01, 1.4420E+01, 1.7174E+01, 2.0735E+01, 2.6381E+01, 3.1959E+01, 3.3691E+01, 3.1933E+01, 2.6295E+01, 2.0551E+01, 1.6872E+01, 1.4005E+01, 1.1797E+01, 9.9605E+00, 8.4376E+00, 7.1582E+00, 6.0408E+00, 5.1470E+00, 4.3476E+00, 3.6516E+00, 3.0404E+00, 2.4991E+00, 2.0178E+00, 1.5875E+00, 1.2065E+00, 8.6981E-01, 5.7090E-01, 3.1174E-01, 8.6615E-02,-1.0611E-01,-2.7020E-01,-4.0798E-01,-5.2420E-01,-6.2164E-01,-7.0014E-01,-7.6514E-01,-8.1616E-01,-8.5783E-01,-8.8987E-01,-9.1467E-01,-9.3309E-01,-9.4505E-01,-9.5353E-01,-9.5628E-01,-9.5692E-01,-9.5673E-01,-9.5382E-01,-9.4210E-01,-9.3682E-01,-9.3215E-01,-9.2239E-01,-9.1062E-01,-8.9979E-01,-8.8954E-01,-8.7715E-01,-8.6582E-01,-8.4632E-01};
static double mSr[mNN]={ 3.0000E+00, 3.6849E+00, 4.2667E+00, 4.8557E+00, 5.4800E+00, 6.1591E+00, 6.9115E+00, 7.7569E+00, 8.7188E+00, 9.8259E+00, 1.1114E+01, 1.2631E+01, 1.4438E+01, 1.6623E+01, 1.9306E+01, 2.2603E+01, 2.7019E+01, 3.2825E+01, 4.2177E+01, 5.1468E+01, 5.4348E+01, 5.1397E+01, 4.1949E+01, 3.2354E+01, 2.6267E+01, 2.1576E+01, 1.8025E+01, 1.5122E+01, 1.2768E+01, 1.0838E+01, 9.1925E+00, 7.9318E+00, 6.8372E+00, 5.9228E+00, 5.1558E+00, 4.5105E+00, 3.9661E+00, 3.5061E+00, 3.1205E+00, 2.7984E+00, 2.5275E+00, 2.3021E+00, 2.1154E+00, 1.9606E+00, 1.8317E+00, 1.7266E+00, 1.6389E+00, 1.5663E+00, 1.5064E+00, 1.4559E+00, 1.4145E+00, 1.3775E+00, 1.3467E+00, 1.3201E+00, 1.2969E+00, 1.2773E+00, 1.2584E+00, 1.2424E+00, 1.2271E+00, 1.2138E+00, 1.2010E+00, 1.1904E+00, 1.1796E+00, 1.1682E+00, 1.1586E+00, 1.1497E+00, 1.1411E+00, 1.1323E+00, 1.1247E+00, 1.1172E+00, 1.1107E+00};
static double mFi[mNN]={ 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 2.9627E-02, 4.5326E-03,-1.2185E-03, 6.2384E-03, 1.9331E-02, 4.8552E-02, 1.0130E-01, 1.9195E-01, 3.2161E-01, 4.0085E-01, 7.4398E-01, 1.0227E+00, 1.3263E+00, 1.6330E+00, 1.9245E+00, 2.1844E+00, 2.4056E+00, 2.5876E+00, 2.7338E+00, 2.8447E+00, 2.9290E+00, 2.9913E+00, 3.0367E+00, 3.0678E+00, 3.0938E+00, 3.1131E+00, 3.1299E+00, 3.1455E+00, 3.1599E+00, 3.1778E+00, 3.1900E+00, 3.2064E+00, 3.2242E+00, 3.2431E+00, 3.2677E+00, 3.2853E+00, 3.3094E+00, 3.3288E+00, 3.3562E+00, 3.3798E+00, 3.4103E+00, 3.4392E+00, 3.4556E+00, 3.4819E+00, 3.5110E+00, 3.5382E+00, 3.5588E+00, 3.5881E+00, 3.6146E+00, 3.6476E+00};
static double mAi[mNN]={ 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 8.3866E+00, 8.2490E+00, 8.0873E+00, 7.8895E+00, 7.6404E+00, 7.3663E+00, 7.0856E+00, 6.8024E+00, 6.5169E+00, 6.1241E+00, 6.0017E+00, 5.7753E+00, 5.5720E+00, 5.3895E+00, 5.2269E+00, 5.0783E+00, 4.9422E+00, 4.8159E+00, 4.6995E+00, 4.5889E+00, 4.4846E+00, 4.3876E+00, 4.2970E+00, 4.2118E+00, 4.1374E+00, 4.0702E+00, 4.0117E+00, 3.9609E+00, 3.9170E+00, 3.8827E+00, 3.8487E+00, 3.8236E+00, 3.8041E+00, 3.7893E+00, 3.7830E+00, 3.7728E+00, 3.7712E+00, 3.7672E+00, 3.7732E+00, 3.7771E+00, 3.7884E+00, 3.8004E+00, 3.8018E+00, 3.8134E+00, 3.8286E+00, 3.8430E+00, 3.8521E+00, 3.8700E+00, 3.8861E+00, 3.9083E+00};
static double mSi[mNN]={ 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 1.4053E+01, 1.3750E+01, 1.3343E+01, 1.2835E+01, 1.2219E+01, 1.1548E+01, 1.0862E+01, 1.0171E+01, 9.4852E+00, 8.6995E+00, 8.2180E+00, 7.6473E+00, 7.1238E+00, 6.6462E+00, 6.2148E+00, 5.8241E+00, 5.4733E+00, 5.1580E+00, 4.8771E+00, 4.6277E+00, 4.4061E+00, 4.2116E+00, 4.0411E+00, 3.8923E+00, 3.7647E+00, 3.6553E+00, 3.5626E+00, 3.4832E+00, 3.4164E+00, 3.3605E+00, 3.3128E+00, 3.2735E+00, 3.2410E+00, 3.2142E+00, 3.1926E+00, 3.1743E+00, 3.1597E+00, 3.1476E+00, 3.1389E+00, 3.1317E+00, 3.1258E+00, 3.1221E+00, 3.1187E+00, 3.1165E+00, 3.1151E+00, 3.1143E+00, 3.1138E+00, 3.1139E+00, 3.1143E+00, 3.1148E+00};



static double Hgam1Fr(double tau) 
{ if(tau< maxTauTab) 
  { double x=1;
    tau-=1;
    if(tau<0) x=1-pow(-tau,0.2);
    else if(tau>0) x=1+pow(tau,0.2);  
    return polint4(x,mNN,mX,mFr); 
  }  else
  {  
     double tau0=maxTauTab;
     return  mFr[mNN-1]
      -log(4*tau) *( log(4*tau)/18  +2./3.) +2*log(tau)   +22/log(4*tau)
    -(-log(4*tau0)*( log(4*tau0)/18 +2./3.) +2*log(tau0)  +22/log(4*tau0));  
  }
}


static double Hgam1Ar(double tau) 
{ if(tau< maxTauTab) 
  { double x=1;
    tau-=1;
    if(tau<0) x=1-pow(-tau,0.2);
    else if(tau>0) x=1+pow(tau,0.2);  
    return polint4(x,mNN,mX,mAr); 
  }  else
  { double tau0=maxTauTab;
    return   mAr[mNN-1]
       -log(4*tau) *( log(4*tau)/18  +2./3.) +2*log(tau)   +21/log(4*tau)
     -(-log(4*tau0)*( log(4*tau0)/18 +2./3.) +2*log(tau0)  +21/log(4*tau0) ); 
  }  
} 

static double Hgam1Sr(double tau) 
{ if(tau< maxTauTab) 
  { double x=1;
    tau-=1;
    if(tau<0) x=1-pow(-tau,0.2);
    else if(tau>0) x=1+pow(tau,0.2);  
    return polint4(x,mNN,mX,mSr); 
  }  else
  { 
     double tau0=maxTauTab;
     return  mSr[mNN-1] + 2/log(4*tau) -2/log(4*tau0); 
  }
}


static double Hgam1Fi(double tau) 
{ 
  if(tau<1) return 0;
  if(tau< maxTauTab) 
  { double x=1;
    tau-=1;
    if(tau<0) x=1-pow(-tau,0.2);
    else if(tau>0) x=1+pow(tau,0.2);  
    return polint4(x,mNN-21,mX+21,mFi+21); 
  }  else 
  {  double tau0=maxTauTab;
     return  mFi[mNN-1] + M_PI/9*(log(4*tau) -log(4*tau0));
  }
}

static double Hgam1Ai(double tau) 
{ 
  if(tau<1) return 0;
  if(tau< maxTauTab) 
  { double x=1;
    tau-=1;
    if(tau<0) x=1-pow(-tau,0.2);
    else if(tau>0) x=1+pow(tau,0.2);  
    return polint4(x,mNN-21,mX+21,mAi+21); 
  }  else 
  {  double tau0=maxTauTab;
     return  mFi[mNN-1] + M_PI/9*(log(4*tau) -log(4*tau0));
  }  
}

static double Hgam1Si(double tau) 
{ 
  if(tau<1) return 0;
  if(tau< maxTauTab) 
  { double x=1;
    tau-=1;
    if(tau<0) x=1-pow(-tau,0.2);
    else if(tau>0) x=1+pow(tau,0.2);  
    return polint4(x,mNN-21,mX+21,mSi+21); 
  }  else  return  mSi[mNN-1];
 
}

double complex Hgam1F(double tau) { return Hgam1Fr(tau)+I*Hgam1Fi(tau); }
double complex Hgam1A(double tau) { return Hgam1Ar(tau)+I*Hgam1Ai(tau); }
double complex Hgam1S(double tau) { return Hgam1Sr(tau)+I*Hgam1Si(tau); }


#ifdef MAIN

int main(void)
{ int i;
  
  double mh=100,mt=171.4, tauT=pow(mh/2/mt,2);
  double GF=1.166637E-5,MW=79.95,SW=0.481,EE=0.31223;
  double tauW=pow(mh/2/MW,2);
  double haa;
  printf("tauT=%e  Hggr(tau)=%E  Hggi(tau)=%E  \n", tauT, HggFr(tauT),HggFi(tauT));
  printf("tauW=%e  HggWr(tauW)=%E  HggWi(tauW)=%E  \n", tauW, HggVr(tauW),HggVi(tauW));
  printf("Hgam1Fr(tauT)=%E\n",Hgam1Fr(tauT));
  printf(" HggFr(tauT)*(1+0.12/M_PI*Hgam1Fr(tauT) )=%E \n", HggFr(tauT)*(1-0.12/M_PI*Hgam1Fr(tauT) ));
  printf("GF=%E(1.166637E-5)\n", sqrt(2.)*pow(EE/(MW*SW),2.)/8.); 
  
//  haa=pow(EE/(MW*SW),2)/8*pow(EE*EE/4/M_PI,2)*pow(mh,3)/(128*pow(M_PI,3))*pow(3*(4./9.)*HggFr(tauT)+HggVr(tauW) ,2);
  haa=pow(EE/(MW*SW),2)/8*pow(EE*EE/4/M_PI,2)*pow(mh,3)/(128*pow(M_PI,3))*pow(HggVr(tauW) ,2);
  printf("haa=%E\n",haa);   
  
//  displayFunc10(mHgam1Fr,-2,4,"mHgam1Fr");    displayFunc10(hdHgam1Fr,-2,4,"hdHgam1Fr");
//  displayFunc10(mHgam1Ar,-0.01,0.01,"mHgam1Ar");    displayFunc10(hdHgam1Ar,-0.01,0.01,"hdHgam1Ar");
//  displayFunc10(mHgam1Sr,-2,4,"mHgam1Sr");    displayFunc10(hdHgam1Sr,-2,4,"hdHgam1Sr");

//  displayFunc10(mHgam1Fi,-2,4,"mHgam1Fi");    displayFunc10(hdHgam1Fi,-2,4,"hdHgam1Fi");
//  displayFunc10(mHgam1Ai,-2,4,"mHgam1Ai");    displayFunc10(hdHgam1Ai,-2,4,"hdHgam1Ai");
//  displayFunc10(mHgam1Si,-2,4,"mHgam1Si");    displayFunc10(hdHgam1Si,-2,4,"hdHgam1Si");
  
//  killPlots();
  
  exit(0);
}

#endif
