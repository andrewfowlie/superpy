
#include "micromegas.h"
#include "micromegas_aux.h"
#include "micromegas_f.h"

double (*sqme22)(int nsub, double GG, double *pvect, int * err_code)=NULL; 
int  nsub22=0;

/*===========================================================*/
static double fixed_Q;
static double GG=1.23;
static double PcmOut, totcoef;
static double pvect[16]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};


static double eps=0.001;

double GGscale=91.187;
/*
double  decayPcm(double am0,  double  am1,  double  am2)
{
  double  summ, diffm, pout;
  summ = am1 + am2;
  diffm = am1 - am2;
  if(am0<summ) return 0;
  return sqrt((am0-summ)*(am0+ summ)*(am0-diffm)*(am0+diffm))/(am0*2);
}
*/          

int  kin22(double PcmIn,double * pmass)
{  
   double sqrtS;
   int i;
   for(i=0;i<16;i++) pvect[i]=0;
   sqrtS=sqrt(pmass[0]*pmass[0]+PcmIn*PcmIn)+sqrt(pmass[1]*pmass[1]+PcmIn*PcmIn);
   PcmOut = decayPcm(sqrtS,pmass[2],pmass[3]);
//printf(" PcmOut =%E (%E %E %E) \n", PcmOut,sqrtS,pmass[2],pmass[3] );   
   if(PcmOut<sqrtS*1.E-4) return 1;
   totcoef = PcmOut /(32.0*M_PI*PcmIn*sqrtS*sqrtS);
   pvect[3] = PcmIn;
   pvect[7] =-PcmIn;
   pvect[0] = sqrt(PcmIn*PcmIn   + pmass[0]*pmass[0]);
   pvect[4] = sqrt(PcmIn*PcmIn   + pmass[1]*pmass[1]);
   pvect[8] = sqrt(PcmOut*PcmOut + pmass[2]*pmass[2]);
   pvect[12]= sqrt(PcmOut*PcmOut + pmass[3]*pmass[3]);

   return 0;
}

double  dSigma_dCos(double  cos_f)
{
   double  r;
   double sin_f=sqrt(fabs((1-cos_f)*(1+cos_f)));
   int err_code=0;
   
   pvect[11]=PcmOut*cos_f;
   pvect[15]=-pvect[11];
   pvect[10]=PcmOut*sin_f;
   pvect[14]=-pvect[10];
   
   
   r = (*sqme22)(nsub22,sqrt(4*M_PI*parton_alpha(GGscale)),pvect,&err_code);
   err_code=0;
   return r * totcoef;
}


double cs22(numout * cc, int nsub, double P, double cos1, double cos2 , int * err) 
{
  int i;
  double pmass[4];
  for(i=1;i<=cc->interface->nvar;i++) if(cc->link[i]) cc->interface->va[i]=*(cc->link[i]);
  GG=sqrt(4*M_PI*parton_alpha(GGscale));
  if( cc->interface->calcFunc()>0 ) {*err=4; return 0;}
  *(cc->interface->gtwidth)=0;
  *(cc->interface->twidth)=0;
  *(cc->interface->gswidth)=1;
  
  for(i=0;i<4;i++) cc->interface->pinf(nsub,1+i,pmass+i,NULL);  
  *err=0;
  sqme22=cc->interface->sqme;
  nsub22=nsub; 
  if(kin22(P,pmass)) return 0; else return 3.8937966E8*simpson(dSigma_dCos,cos1,cos2,0.3*eps);
}

/*===================  Collider production ==========*/

static double sMin,sMax,pcmOut;
static double pmass[4];
static double x0;
static int pc1_,pc2_;
static int ppFlag;

static numout * colliderProduction(char * name1,char *name2)
{ 
  char libname[100], process[100], lName1[20], lName2[20];
  numout *cc;
  int i,first;
    
  pname2lib(name1,lName1);
  pname2lib(name2,lName2);
  if(strcmp(lName1,lName2)>0)sprintf(libname,"PP_%s%s",lName1,lName2);
  else                       sprintf(libname,"PP_%s%s",lName2,lName1); 
  sprintf(process,"proton,proton->%s,%s{",name1,name2);
  
  for(i=0,first=1;i<nModelParticles;i++)
  { switch(abs(ModelPrtcls[i].NPDG))
    { case 1: case 2: case 3: case 4: case 5: case 81: case 83:
       if(!first) strcat(process,","); else  first=0;
       sprintf(process+strlen(process),"%s,%s", 
                     ModelPrtcls[i].name,ModelPrtcls[i].aname);  
       break;
       case 21:
       if(!first) strcat(process,",");else  first=0;
       sprintf(process+strlen(process),"%s",
                             ModelPrtcls[i].name);
    }                              
  }
  
  cc=getMEcode(0,ForceUG,process,NULL,"",libname);
  

  return cc;
}


static double  cos_integrand(double xcos)
{ int err;
  double xsin=sqrt(1-xcos*xcos);
  double GG;
  pvect[9]=pcmOut*xcos;
  pvect[10]=pcmOut*xsin;
  pvect[13]=-pvect[9];
  pvect[14]=-pvect[10];
/*  
  s=(pvect[0]+pvect[4]); s=s*s;
  t=pmass[0]*pmass[0]+pmass[2]*pmass[2]-2*(pvect[0]*pvect[8] -pvect[1]*pvect[9]);
  u=pmass[0]*pmass[0]+pmass[3]*pmass[3]-2*(pvect[0]*pvect[12]-pvect[1]*pvect[13]);
  Q=sqrt(2*s*fabs(t*u)/(s*s+t*t+u*u));
*/

  GG=sqrt(4*M_PI*parton_alpha(fixed_Q));
  return  sqme22(nsub22,GG,pvect,&err)*convStrFun2(x0,fixed_Q,pc1_,pc2_,ppFlag);  
}


static double  s_integrand(double y)
{  double r;
   double pp=6;
   double s,pcmIn;
   s=sMin+pow(y,pp)*(sMax-sMin);
   
   pcmIn=decayPcm(sqrt(s),pmass[0], pmass[1]);
   if(pcmIn==0) return 0;
   pvect[0]=sqrt(pmass[0]*pmass[0]+pcmIn*pcmIn);
   pvect[1]=pcmIn; pvect[2]=0; pvect[3]=0;
   pvect[4]=sqrt(pmass[1]*pmass[1]+pcmIn*pcmIn);
   pvect[5]=-pcmIn; pvect[6]=0; pvect[7]=0;
   pcmOut=decayPcm(sqrt(s),pmass[2], pmass[3]);
   pvect[8]=sqrt(pmass[2]*pmass[2]+pcmOut*pcmOut);
   pvect[11]=0;
   pvect[12]=sqrt(pmass[3]*pmass[3]+pcmOut*pcmOut);
   pvect[15]=0;
   x0=s/sMax;
 
   r=  3.8937966E8*pcmOut/(32*M_PI*pcmIn*s)*simpson(cos_integrand,-1.,1.,1.E-3);
   r*=pp*pow(y,pp-1)*(sMax-sMin);
   return r; 
}

double hCollider(double Pcm, int pp,int sc, char * name1,char *name2)
{ 
  double  sigma_tot=0, Qstat;
  int i;
  numout *cc;
  int n1,n2;
  
  n1=pTabPos(name1);  if(n1==0) { printf("%s - no such particle\n",name1); return 0;}
  n2=pTabPos(name2);  if(n2==0) { printf("%s - no such particle\n",name2); return 0;}
      
 
  sMax=4*Pcm*Pcm; 
  sMin=pMass(name1)+pMass(name2); sMin*=sMin; sMin+=1;
  ppFlag=pp;   
  cc=colliderProduction( name1,name2);
  if(!cc) return 0;


  for(i=1;i<=cc->interface->nvar;i++) 
  { if(cc->link[i]) cc->interface->va[i]=*(cc->link[i]);
  } 
  
  fixed_Q=sqrt(sMin)/2;
  if(sc) 
  {  if( ModelPrtcls[abs(n1)-1].cdim==8 && ModelPrtcls[abs(n2)-1].cdim==8 )
     fixed_Q*=0.1; else   fixed_Q*=0.2;
  }
  if(Qaddress)
  { Qstat=*Qaddress;
    *Qaddress=fixed_Q;
  }  
  if( cc->interface->calcFunc()>0 )  return 0;
  *(cc->interface->gtwidth)=0;
  *(cc->interface->twidth)=0;
  *(cc->interface->gswidth)=0;
  sqme22=cc->interface->sqme;
      
  sigma_tot=0;   
  for(nsub22=1;nsub22<=cc->interface->nprc; nsub22++) 
  { int pc[4];
    double tmp;

    for(i=0;i<4;i++) cc->interface->pinf(nsub22,i+1,pmass+i,pc+i);

    if(pc[0]<=pc[1])
    { pc1_=pc[0];
      pc2_=pc[1];
      tmp=simpson(s_integrand,0.,1.,1.E-2)/sMax;     
      sigma_tot+=tmp;
    }  
  }
  if(Qaddress){ *Qaddress=Qstat;}
            
  return sigma_tot;
}

#ifdef plazmaWidth
static numout* plazmaWidth_cc;
static double plazmaWidth_T;
static double plazmaWidth_m[4];
static double plazmaWidth_integrand(double Pcm)
{ int err;
  double E1,E2,sqrt_s; 
  if(Pcm==0) return 0;  
  E1=sqrt(Pcm*Pcm+plazmaWidth_m[0]*plazmaWidth_m[0]);
  E2=sqrt(Pcm*Pcm+plazmaWidth_m[1]*plazmaWidth_m[1]);
  sqrt_s=E1+E2;
  if(sqrt_s<=plazmaWidth_m[2]+plazmaWidth_m[3]) return 0;
  
  return 4*bessk1(sqrt_s/plazmaWidth_T)*cs22(plazmaWidth_cc,1,Pcm, -1., 1. , &err)*pow(sqrt_s*Pcm,3.)/E1/E2; 
}

double plazmaWidth(char *process,double T)
{  char libName[40];
   plazmaWidth_T=T;
   process2Lib(process,libName);
   process2Mass(process,plazmaWidth_m);
   plazmaWidth_cc=getMEcode(0,ForceUG,process,NULL,NULL,libName); 
   return simpson(plazmaWidth_integrand,0., 5*T,1.E-3)/(4*M_PI*M_PI*plazmaWidth_m[1]*plazmaWidth_m[1]*bessk2(plazmaWidth_m[1]/T))
   /3.8937966E8;    
}
#endif
/*============ Fortran ==========*/

double cs22_(int*ccf,int*nsub,double*P,double*cos1,double*cos2,int*err)
{ numout*cc;
  memcpy(&cc,ccf,sizeof(cc));
  return cs22(cc,*nsub,*P,*cos1,*cos2 ,err);
} 
void  sethelicities_(double *h1,double *h2) { Helicity[0]=*h1; Helicity[1]=*h2;} 

