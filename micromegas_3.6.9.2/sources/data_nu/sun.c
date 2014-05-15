#include"../sources/micromegas.h"
#include"../sources/micromegas_aux.h"
#include"lib/pmodel.h"

double S(double x, double v)
{double A,Ab;
double vesc=1156;  /* escape velocity for the Sun */
A=3*x/2./(x-1)/(x-1)*vesc*vesc/v/v;
Ab=pow(A,1.5);
return pow((Ab/(1+Ab)),2./3.);
}



double capture(double Mwimp, double csSIp,double csSIn,double csSDp,double csSDn,double v )
{ 
   double mp=0.938;
/* sum over nuclear species */
/* H,He,C,N,O,Ne,Mg,Ni,S,Fe  */
double A[10]={1,4,12,14,16,20,24,28,32,56};
double Z[10]={1,2,6,7,8,10,12,14,16,26};
double fs[10]={0.772,0.209,3.87e-3,9.4e-4,8.55e-3,1.51e-3,7.39e-4,8.13e-4,4.65e-4,1.46e-3};
double fsp[10]={0.670,0.311,2.37e-3,1.88e-3,8.78e-3,1.93e-3,7.33e-4,7.98e-4,5.50e-4,1.42e-3};

double ps[10]={3.16,3.4,3.23,3.23,3.23,3.23,3.23,3.23,3.23,3.23};
double psp[10]={3.15,3.4,2.85,3.83,3.25,3.22,3.22,3.22,3.22,3.22};
double Finf[10]={1,0.986,0.788,0.613,0.613,0.613,0.281,0.281,0.101,0.00677};
double mc[10]={1,18.2,61.6,75.2,75.2,75.2,71.7,71.7,57.0,29.3};
double ai[10]={1,1.58,2.69,2.69,2.69,2.69,2.97,2.97,3.1,3.36};
double si[10],mu[10];
double sum,cap,fp,fn,form,sumsi,ap,an;
int i;
fp=sqrt(csSIp)*(Mwimp+mp)/2/Mwimp/mp;
fn=sqrt(csSIn)*(Mwimp+mp)/2/Mwimp/mp;

printf("fp=%.3e fn=%.3e\n",fp,fn);
sum=0;
sumsi=0.;
for(i=0;i<10;i++)
{mu[i]=Mwimp*Mwimp*mp*mp*A[i]*A[i]/(Mwimp+A[i]*mp)/(Mwimp+A[i]*mp);
si[i]=pow((Z[i]*fp+(A[i]-Z[i])*fn),2)*1.e4*4.*mu[i];
form=(Finf[i]+(1-Finf[i])*exp(-pow((log(Mwimp)/log(mc[i])),ai[i])));

sum+=(Finf[i]+(1-Finf[i])*exp(-pow((log(Mwimp)/log(mc[i])),ai[i])))*fs[i]*ps[i]/A[i]/mp*
S(Mwimp/A[i]/mp,v)*si[i];
sumsi+=(Finf[i]+(1-Finf[i])*exp(-pow((log(Mwimp)/log(mc[i])),ai[i])))*fsp[i]*psp[i]/A[i]/mp*
S(Mwimp/A[i]/mp,v)*si[i];
/*printf("i=%d Z=%.2e A=%.2e si_he=%.3e mu=%.2e fo=%.2e sum=%.3e sump=%.3e\n",i,Z[i],A[i],si[i],mu[i],form,sum,sump);
*/}


return cap=4.8e24/Mwimp*sumsi+1.3e25*csSDp*1.e4*S(Mwimp/mp,v)/Mcdm;

}

/*
int main(int narg, char** args)
{ double csSIp,csSIn,csSDp,csSDn;
  double vRotation=220, vSun=225.2,vEsc=600.;

  if(narg < 6)
  {
    printf("Parameters expected: Mwinp[GeV], csSIp, csSIn, csSDp,csSDn[pb]\n"); 
    return 1;
  }   

  if(sscanf(args[1],"%lf",&Mcdm)==1 &&
     sscanf(args[2],"%lf",&csSIp)==1 &&
     sscanf(args[3],"%lf",&csSIn)==1 &&
     sscanf(args[4],"%lf",&csSDp)==1 &&
     sscanf(args[5],"%lf",&csSDn)==1) 
    printf("Model parameters:\n Mwimp=%.2E csSIp=%.2E, csSIn=%.2E,csSDp=%.2E,csSDn=%.2E\n",
                          Mcdm, csSIp, csSIn,csSDp,csSDn);
  else { printf("Wrong parameters\n"); return 2;}  
 
 
  SetfMaxwell(vRotation,vEsc); 
  capture(Mcdm, csSIp,csSIn,csSDp,csSDn,vRotation );
 
  return 0;

}
*/
  
  
  
  
  
  
  
  
  
  

 
