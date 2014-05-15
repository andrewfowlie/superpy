#include "SLHAplus.h"
 
double MbPole=4.967923;  
static double MtPole=174.1; 
static double qMin=3.125347E-1, qLim=1;
static double m_fact(int nf, double alpha1, double alpha2);
static int nfMax=6;
/*
   this file contains routines for running QCD coupling, 
   running quark masses,  and  effective quarks masses 
   which  generate correct widths of scalar particles.
*/


double Mbp(void) { return MbPole;}


static double alpha3(double Q, double lambda, int nf)
{
  double b0,b1,b2,lg,lg2,x;

  b0 =   11 -    (2./3.)*nf;
  b1 =   51 -   (19./3.)*nf;
  b2 = 2857 - (5033./9.)*nf + (325./27.) *nf*nf;
  lg = 2*log(Q/lambda);
  lg2= log(lg);
  x  = 2*b1/(b0*b0*lg);

  return 4*M_PI*(1-x*(lg2 -x*((lg2-0.5)*(lg2-0.5)+b2*b0/(8*b1*b1)-1.25) ))/(b0*lg);
}


static double findLambda(int nf,double alpha, double M)
{ double l1=0.1, l2=0.3;
  double l,a,a1,a2;

  l2=M*exp(-2*M_PI/alpha/(11-2*nf/3.));
  while((a2=alpha3(M,l2,nf)-alpha) < 0) l2*=1.2;
  l1=l2*0.8;
  while((a1=alpha3(M,l1,nf)-alpha) > 0)  l1*=0.8;

  do{ l=(l1*a2-l2*a1)/(a2-a1);  
      a=alpha3(M,l,nf)-alpha;
      if(a<0) { a1=a;l1=l;} else {a2=a;l2=l;}
    } while (fabs(a) > 0.00001*alpha);
  return l;
}      

static double qMass[7]={0,0,0,0, 1.200000E+00, 4.230000E+00, 1.619152E+02};
static double lambda[7]={0,0,0,3.125347E-01, 2.763267E-01, 1.991586E-01, 8.449407E-02};


double poleQmass(double mm, double alpha, int nf)
{
  double 
         zeta2=1.6439,
         c1=4./(3.),
         c2=13.4434-1.0414*(nf-1),
         c3=190.595-(nf-1)*(26.655-(nf-1)*0.6527);
  double a=alpha/M_PI;

       if(nf==5) c2+=zeta2*(1.61/4.62); /* (Mc+Ms)/Mb contribution */
  else if(nf==6) c2+=zeta2*(6.23/175); /* (Mc+Ms+Mb)/Mt contribution */

  return mm*(1+a*(c1 +a*(c2 +a*c3)));   

}

static int notInitialized=1;

double  initQCD(double MZalphaS, double McMc, double MbMb, double MtP)
{ 
  double Mq,Mq_;

  lambda[5]= findLambda(5,MZalphaS, 91.187);

  MtPole=MtP;

  if(MtP< poleQmass(91.187,MZalphaS, 6)) return -1;  
  
  for(Mq=MtP,Mq_=0;fabs(Mq_-Mq)>0.00001*Mq;)
  { double alpha=alpha3(Mq, lambda[5], 5);
    Mq_=Mq; 
    Mq*=MtP/poleQmass(Mq, alpha, 6);
  }

  qMass[6]= Mq;

  lambda[6]= findLambda(6,alpha3(qMass[6],lambda[5], 5),qMass[6]);
  notInitialized=0;
  nfMax=6;
  
  qMass[5]=0; qMass[4]=0;
  if(MbMb<=lambda[5]) { qMin=lambda[5]; return qMin;}
  
  qMass[5]=MbMb;
  MbPole=poleQmass(MbMb, alpha3(MbMb,lambda[5] ,5),5);

  lambda[4]= findLambda(4,alpha3(qMass[5],lambda[5],5),qMass[5]);

  if(McMc<=lambda[4]) {qMin=lambda[4]; return qMin;}
  qMass[4]=McMc;
  lambda[3]=findLambda(3,alpha3(qMass[4],lambda[4],4),qMass[4]);
                       qMin=lambda[3]; return qMin;
}

double  initQCD5(double MZalphaS, double McMc, double MbMb, double MtP)
{
   double lmbd= initQCD(MZalphaS,McMc, MbMb,MtP);
   lambda[6]=lambda[5]; 
   nfMax=5;
   return lmbd;
}


static int  NF(double Q)
{ 
         if(Q<qMass[4]) return 3;
   else  if(Q<qMass[5]) return 4; 
   else  if(Q<qMass[6]||nfMax==5) return 5; 
   else                 return 6;
}

double alphaQCD(double Q) 
{ 
  if(notInitialized) initQCD(0.1184,1.2,4.23,173.07);
  if(Q<qLim) Q=qLim;
  if(Q<qMin) return 1; return alpha3(Q,lambda[NF(Q)],NF(Q));
}



static double m_fact(int nf, double alpha1, double alpha2)
{  double k,c1,c2;

   switch(nf)
   { case 3: k=4./9;   c1=0.895; c2=1.371; break;
     case 4: k=12./25; c1=1.014; c2=1.389; break;
     case 5: k=12./23; c1=1.175; c2=1.501; break;
     case 6: k=4./7;   c1=1.398; c2=1.793; break;
   }

  { double xx=
           pow(2*alpha2/k,k)*(1+alpha2*(c1+alpha2*c2))/
           (pow(2*alpha1/k,k)*(1+alpha1*(c1+alpha1*c2)));
    return xx;
  }          
}
  

static double  runningMass(double M0, double Q0,double Q)
{  double alpha0, alpha1;
   int n=0;   
   if(Q<qLim) Q=qLim;
   for(n=3; n<6;n++) if(Q0<qMass[n+1]) break;     

   alpha0=alphaQCD(Q0)/M_PI;   

   if(Q<Q0)
     for(;Q<Q0; n--,alpha0=alpha1 )
     { 
        if(Q<qMass[n]) Q0=qMass[n]; else Q0=Q;
        alpha1=alphaQCD(Q0)/M_PI;
        M0*=m_fact(n,alpha0,alpha1);     
     }
   else 
     for(;Q>Q0;n++,alpha0=alpha1)                    
     {  
        if(n<6 && Q>qMass[n+1]) Q0=qMass[n+1]; else Q0=Q;
        alpha1=alphaQCD(Q0)/M_PI;
        M0*=m_fact(n,alpha0,alpha1);
     }    
   return M0;
}

double MqRun(double mass2GeV, double Q)
{ return runningMass(mass2GeV, 2., Q);}



double McRun(double Q) 
{  if(Q<=qMass[4])  return qMass[4];
   if(qMass[4]==0) return 0;
   return runningMass(qMass[4], qMass[4], Q); 
}

double MbRun(double Q) 
{  if(Q<qMass[5]) return qMass[5];
   if(qMass[5]==0) return 0; 
   return runningMass(qMass[5], qMass[5], Q); 
}

double MtRun(double Q) 
{
   return runningMass(qMass[6], qMass[6], Q); 
} 

static double DeltaQCD(double Q)
{
 double  a=alphaQCD(Q)/M_PI;
 int nf=NF(Q);
 double res, res1;

 return  a*( 5.67 + a*( (35.94-1.36*nf) + a*(164.14-nf*(25.77-nf*0.259) )));
}

double MqEff(double mass2GeV, double Q) { return MqRun(mass2GeV,Q)*sqrt(1+DeltaQCD(Q));}

double McEff(double Q) { return McRun(Q)*sqrt(1+DeltaQCD(Q)); }
double MbEff(double Q) { double m=MbRun(Q)*sqrt(1+DeltaQCD(Q));  if(m>MbPole) return MbPole; else return m;}

double MtEff(double Q) 
{ 
  double meff,alpha;
  if(Q<=2*MtPole) return  MtPole;
  meff=MtRun(Q)*sqrt(1+DeltaQCD(Q));
  if(meff>MtPole)   return  MtPole;
  alpha=pow(2*MtPole/Q,15.);
  return alpha*MtPole+(1-alpha)*meff;
}

 
double nfQCD(double Q) {return NF(Q);}


#ifdef TEST
int main(int n, char **args)
{
  double alphaSMZ,mbp,mbmb,mtp;
  double Q; 
  int nf;
  double McMc=1.4;
  sscanf(args[1],"%lf",&alphaSMZ);
  sscanf(args[2],"%lf",&McMc);
  sscanf(args[3],"%lf",&mbmb);
  sscanf(args[4],"%lf",&mtp);
  sscanf(args[5],"%lf",&Q);

  initQCD(alphaSMZ,McMc,mbmb,mtp);

printf("qMass : %E %E %E %E\n",  qMass[3], qMass[4], qMass[5], qMass[6]); 
printf("lambda: %E %E %E %E\n", lambda[3],lambda[4],lambda[5],lambda[6]);
printf("qMin=%E\n",qMin);

  printf("MbPole=%f\n", MbPole);
  printf("MtMt=%f\n", qMass[6]);

  printf("qmass[6]=%f\n",qMass[6]);

  printf("alphaS(%f)=%f\n",Q,alphaQCD(Q));
  printf("MbRun=%f  MbEff=%f \n",MbRun(Q),MbEff(Q)  );
  printf("MtRun=%f  MtEff=%f \n",MtRun(Q),MtEff(Q)  );
}
#endif
