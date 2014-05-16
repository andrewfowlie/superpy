#include"../sources/micromegas.h"
#include"../sources/micromegas_aux.h"
#include"lib/pmodel.h"

double SE(double x, double v)
{double A,Ab;
double vesc=13.2;  /* escape velocity for the Earth */
A=3*x/2./(x-1)/(x-1)*vesc*vesc/v/v;
Ab=pow(A,1.5);
return pow((Ab/(1+Ab)),2./3.);
}



double captureE(double Mwimp, double csSIp,double csSIn,double csSDp,double csSDn,double v )
{ 
   double mp=0.938;
/* sum over nuclear species */
/* H,He,C,N,O,Ne,Mg,Ni,S,Fe  */
double A[9]={16,28,24,56,40,30,23,32,59};
double Z[9]={8,14,12,26,20,15,11,16,28};
double fs[9]={0.3,0.15,0.14,0.3,0.015,0.011,0.004,0.05,0.03};
double ps[9]={1.2,1.2,1.2,1.6,1.2,1.2,1.2,1.6,1.6};
double si[9],mu[9],FE[9];
double sum,cap,fp,fn,form,sumsi,ap,an,Ax,x;
double vesc=13.2;  /* escape velocity for the Earth */
double rho=0.3;
int i;
fp=sqrt(csSIp)*(Mwimp+mp)/2/Mwimp/mp;
fn=sqrt(csSIn)*(Mwimp+mp)/2/Mwimp/mp;

printf("fp=%.3e fn=%.3e\n",fp,fn);
sum=0;

for(i=0;i<10;i++)
{mu[i]=Mwimp*Mwimp*mp*mp*A[i]*A[i]/(Mwimp+A[i]*mp)/(Mwimp+A[i]*mp);
si[i]=pow((Z[i]*fp+(A[i]-Z[i])*fn),2)*1.e4*4.*mu[i];
 FE[i]=1;
 if(i==3)
 { x=Mwimp/A[i]/mp;
 Ax=3*x/2./(x-1)/(x-1)*vesc*vesc/v/v;
 FE[i]=1-0.26*(Ax)/(1+Ax);}
sum+=FE[i]*fs[i]*ps[i]/A[i]/mp*
SE(Mwimp/A[i]/mp,v)*si[i];

printf("i=%d Z=%.2e A=%.2e si_he=%.3e mu=%.2e fo=%.2e sum=%.3e sump=%.3e\n",i,Z[i],A[i],si[i],FE[i],sum);
}

printf("si=%.3e \n",4.8e15/Mwimp*sum);
printf("S=%.3e\n",SE(Mwimp/mp,v));

return cap=4.8e15/Mwimp*sum/(v/270.)*rho/0.3;

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
  
  
  
  
  
  
  
  
  
  

 
