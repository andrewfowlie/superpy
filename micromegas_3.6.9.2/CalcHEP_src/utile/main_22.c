#include<math.h>
#include<stdio.h>
#include<unistd.h>
#include<sys/stat.h>
#include<sys/types.h>

#include <dlfcn.h>
#include <sys/wait.h> 
             

#include"nType.h"
#include"num_in.h"
#include"num_out.h"
#include"VandP.h"
#include"dynamic_cs.h"
#include"rootDir.h" 


int main(void)
{ int err;

  REAL pvec[16];
  REAL *p1=pvec, *p2=pvec+4, *p3=pvec+8, *p4=pvec+12;

  REAL m1,m2,m3,m4,m[4];

  double  Pin,Pout;
  double  totcoef, sqrt_S, S, lambda12, lambda34,ms,md;
  double  cos_fi, cos_step;
  double  sigmaTot=0;
  int i;
  char mess[20];
  double GG;
  numout*cc;
  txtList decays;

// Specify model for work: Model directory and model number.
  setModel("models" ,2 ); 

// Calculation of public constraints  
  err=calcMainFunc();
  if(err) { printf("Can not calculate constrained parameter %s\n",varNames[err]);return err;}
  
  printf("GF=%E\n",findValW("GF"));
  slhaDecayPrint("W+", stdout);

// SQME code generation 
   cc= newProcess("e,E->m,M");
   err=passParameters(cc);
   if(err) { printf("Can not calculate constrained parameter %s\n",cc->interface->varName[err]); return err;}

   Pin=100;
 
/* find masses */ 
  procInfo2(cc,1,NULL, m);
  m1=m[0];m2=m[1];m3=m[2];m4=m[3];

/* 2-2 kinematics */

  sqrt_S=sqrt(m1*m1+Pin*Pin) +sqrt(m2*m2+Pin*Pin);
  S=sqrt_S*sqrt_S;  
  lambda12=2*sqrt_S*Pin;
  ms = m3+m4; if (ms >= sqrt_S) return 1;
  md = m3 -m4;
  lambda34 = sqrt((S - ms*ms) * (S - md*md));
  totcoef = 3.8937966E8 * lambda34 /(32.0 * M_PI * lambda12 * S);        
  Pout=lambda34/(2*sqrt_S);

/* fill momenta of particles */
  
  for(i=0;i<16;i++) pvec[i]=0;
  
  p1[0]=  sqrt(Pin*Pin + m1*m1);

  
  p1[3]=  Pin;
  p2[0]=  sqrt(Pin*Pin + m2*m2);
  p2[3]= -Pin;
  p3[0]=  sqrt(Pout*Pout + m3*m3);
  p4[0]=  sqrt(Pout*Pout + m4*m4);
  
/* Assign QCD coupling */
  GG=1.238;
  
  
  cos_step=0.1; /* step to calculate dsigma/dcos */
  sigmaTot=0.;  /* total cross section */
/* cycle */
  for(cos_fi=1-cos_step/2; cos_fi>-1; cos_fi -= cos_step)
  { double  DsigmaDcos;
    double  sin_fi=sqrt((1-cos_fi)*(1+cos_fi));
    int err=0;
    p3[3]=Pout*cos_fi; p4[3]=-p3[3];
    p3[2]=Pout*sin_fi; p4[2]=-p3[2];
    DsigmaDcos=totcoef*cc->interface->sqme(1,GG,pvec,&err); 
    sigmaTot += DsigmaDcos*cos_step;
  }
  printf("sigmaTot(Pcm=%E)= %E\n",Pin,sigmaTot);

  return 0;
}
