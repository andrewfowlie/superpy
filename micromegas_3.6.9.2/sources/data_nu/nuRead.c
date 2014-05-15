/*
Earth_nueb_EVOL.dat
Earth_nue_EVOL.dat
Earth_numub_EVOL.dat
Earth_numu_EVOL.dat
Earth_nutaub_EVOL.dat
Earth_nutau_EVOL.dat
explanations.r3-1.pdf
*/

/*
Sun_nueb_EVOL.dat
Sun_nue_EVOL.dat
Sun_numub_EVOL.dat
Sun_numu_EVOL.dat
Sun_nutaub_EVOL.dat
Sun_nutau_EVOL.dat
*/

#include<stdio.h>
#include<stdlib.h>
#include"../sources/micromegas.h"
#include"../sources/micromegas_aux.h"

double massTabSun[12]={10, 30, 50, 70, 90, 100, 200, 300, 500, 700, 900, 1000};

double xTab[100]={0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.2,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,0.3,0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.39,0.4,0.41,0.42,0.43,0.44,0.45,0.46,0.47,0.48,0.49,0.5,0.51,0.52,0.53,0.54,0.55,0.56,0.57,0.58,0.59,0.6,0.61,0.62,0.63,0.64,0.65,0.66,0.67,0.68,0.69,0.7,0.71,0.72,0.73,0.74,0.75,0.76,0.77,0.78,0.79,0.8,0.81,0.82,0.83,0.84,0.85,0.86,0.87,0.88,0.89,0.9,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99,1.};
double nuTabSun[12][2][8][100];
char * chan[8]={"nu","b","tau","c","q","t","Z","W"};
double nuTabEarth[14][2][8][100];
double massTabEarth[14]= {10, 30, 50, 70, 90, 100, 150, 200, 250, 300, 500, 700, 900, 1000};

static int OK=0;

int   nuRead(void)
{ int i,j,k,l,nLine=2;
  double x,mass;

  for(l=0;l<2;l++)
  {
     FILE*f;
     if(l==0) f=fopen("Sun_numu_EVOL.dat","r"); else f=fopen("Sun_numub_EVOL.dat","r");
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

  for(l=0;l<2;l++)
  {
     FILE*f;
     if(l==0) f=fopen("Earth_numu_EVOL.dat","r"); else f=fopen("Earth_numub_EVOL.dat","r");
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

  OK=1;
  return 0;  
}


static void mInterpSun(double Nmass,  int  CHin, int  CHout, double*tab)
{  
   int l,i0;
   double c0,c1;
   double *p0,*p1;
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
  
  
                                                   
int basicNuSpectraSun(int pdgN, int outN, double * tab)
{ 
  int inP,i,j;
  int N=abs(pdgN);
  double tab100[100];
  double c=0.5;
  
  for(i=0;i<NZ;i++) tab[i]=0;
  
  if(N==12 ||N==14 ||N==16) {  if( pdgN*outN<0) return 0; else c=1;}  
  
  switch(N)
  {
    case 1: case 2: case 3: inP=4; break;
    case 4:                 inP=3; break;
    case 5:                 inP=1; break; 
    case 6:                 inP=5; break;
 
    case 15:                inP=2; break;  /*l  */ 
    case 23:                inP=6; break;  /*z  */ 
    case 24:                inP=7; break;  /*w  */
    case 12: 
    case 14: 
    case 16:                inP=0; break;  /*nu */
    default:                       return 1;
  }  
  if(!OK) nuRead();

  mInterpSun(Mcdm,inP,-(outN-1)/2,tab100);

  for(i=0;i<NZ;i++)
  { double x=exp(Zi(i));
    tab[i]+=c*x*polint3(x,100, xTab ,tab100);  
  }

  return 0;
}


int basicNuSpectraEarth(int pdgN, int outN, double * tab)
{ 
  int inP,i,j;
  int N=abs(pdgN);
  double tab100[100];
  double c=0.5;
  
  for(i=0;i<NZ;i++) tab[i]=0;
  
  if(N==12 ||N==14 ||N==16) {  if( pdgN*outN<0) return 0; else c=1;}
  
  switch(N)
  {
    case 1: case 2: case 3: inP=4; break;
    case 4:                 inP=3; break;
    case 5:                 inP=1; break; 
    case 6:                 inP=5; break;
 
    case 15:                inP=2; break;  /*l  */ 
    case 23:                inP=6; break;  /*z  */ 
    case 24:                inP=7; break;  /*w  */
    case 12:
    case 14:
    case 16:             tab[0]=2/(Zi(0)-Zi(1));
                            return 0;     /*nu */
    default:                       return 1;
  }  
  if(!OK) nuRead(); 
  mInterpEarth(Mcdm,inP,-(outN-1)/2,tab100);
  for(i=0;i<NZ;i++)
  { double x=exp(Zi(i));
    tab[i]+=c*x*polint3(x,100, xTab ,tab100);  
  }

  return 0;
}


static int PrintOn=0;
static long stdPDG(long N)
{ 
   switch (N)
   { case  81: return  1;
     case -81: return -1;
     case  83: return  3;
     case -83: return -3;
     default : return N; 
   }    
}       
        
/* v*cs22 at v=0 */
static double  vcs22(numout * cc,int nsub,int * err)
{
   int i;
   double pcm,r;
   double pmass[4], pvect[16];
   
   for(i=1;i<=cc->interface->nvar;i++) if(cc->link[i]) 
   cc->interface->va[i]=*(cc->link[i]);

   if( cc->interface->calcFunc()>0 ) {*err=4; return 0;}
   *(cc->interface->gtwidth)=0;
   *(cc->interface->twidth)=0;
   *(cc->interface->gswidth)=0;
   for(i=0;i<4;i++) cc->interface->pinf(nsub,1+i,pmass+i,NULL);
   *err=0;
   if(pmass[0]+pmass[1] <= pmass[2]+pmass[4]) return 0;
   for(i=0;i<16;i++) pvect[i]=0;

   pcm= decayPcm(pmass[0]+pmass[1],pmass[2],pmass[3]);
   for(i=0;i<2; i++) pvect[4*i]=pmass[i];
   for(i=2;i<4; i++) pvect[4*i]=sqrt(pmass[i]*pmass[i] +pcm*pcm);
   pvect[8+3]=pcm;
   pvect[12+3]=-pcm;
   r=cc->interface->sqme(nsub,pvect,err);

   return 3.8937966E8*r*pcm/(16*M_PI*pmass[0]*pmass[1]*(pmass[0]+pmass[1]));
}


static double calcSpectrum0(char *name1,char*name2, int key, double *SpectranuSun, double *SpectraNuSun, double *SpectranuEarth, double *SpectraNuEarth )
{
  int i,k;
  double vcsSum=0; 
  int ntot,err;
  double * v_cs;
  
  char name1L[10],name2L[10], lib[20],process[20];
  numout * libPtr;
  
  for(i=0;i<NZ;i++) SpectranuSun[i]=SpectraNuSun[i]=SpectranuEarth[i]=SpectraNuEarth[i]=0;  

  pname2lib(name1,name1L);
  pname2lib(name2,name2L);
  sprintf(lib,"omg_%s%s",name1L,name2L);
  sprintf(process,"%s,%s->2*x",name1,name2);
  libPtr=getMEcode(0,ForceUG,process,NULL,txtListOddParticles(),lib); 
  if(!libPtr) return 0;

  procInfo1(libPtr,&ntot,NULL,NULL); 

  
  v_cs=malloc(sizeof(double)*ntot);
  (*libPtr->interface->twidth)=0;
  
  for(k=0;k<ntot;k++)
  { double m[4];
    char *N[4];
    procInfo2(libPtr,k+1,N,m);
    if((m[2]+m[3])/(m[0]+m[1])<1)
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
   
  for(k=0;k<ntot ;k++) if(v_cs[k]>=0)
  { char * N[4];
    double m[4];
    int l, charge3[2],spin2[2],cdim[2],pdg[2];
    int PlusAok=0;

    procInfo2(libPtr,k+1,N,m);
    for(l=0;l<2;l++)  pdg[l]=stdPDG(qNumbers(N[2+l],spin2+l,charge3+l,NULL));

            
    if(v_cs[k]>1.E-3*vcsSum) 
    {  double tab2[NZ]; 
       if(PrintOn )
       { char txt[100];
         sprintf(txt,"%s,%s -> %s %s", N[0],N[1],N[2],N[3]);
         printf("  %-20.20s  %.2E\n",txt,v_cs[k]*2.9979E-26);
       }
       for(l=0;l<2;l++)   
       {      
         basicNuSpectraSun(pdg[l],1,tab2);
          for(i=0;i<NZ;i++) SpectranuSun[i]+=tab2[i]*v_cs[k];
         basicNuSpectraSun(pdg[l],-1,tab2);
          for(i=0;i<NZ;i++) SpectraNuSun[i]+=tab2[i]*v_cs[k];  
         basicNuSpectraEarth(pdg[l],1,tab2);
          for(i=0;i<NZ;i++) SpectranuEarth[i]+=tab2[i]*v_cs[k];
         basicNuSpectraEarth(pdg[l],-1,tab2);
          for(i=0;i<NZ;i++) SpectraNuEarth[i]+=tab2[i]*v_cs[k];  
       }
    } 
  }
  free(v_cs);
  return  vcsSum;
}



double calcNuSpectrum(int key, double*SnmSun, double*SNmSun, double*SnmEarth, double*SNmEarth, int *errcode)
{ int n,i,l,err;
  double vcs;
  numout* cc;
  char  lop[20];
  double *Spectra,*Spectra_;
  double  buffnSun[NZ],buffNSun[NZ],buffnEarth[NZ],buffNEarth[NZ];
  char * name, *aname;
  
  if(key&4) { PrintOn=1; printf("    Channel          vcs[cs^3/s]\n");} else PrintOn=0;  
  vcs=0;
  if(errcode) *errcode=0;  

  err=sortOddParticles(lop); 
  if(err) { printf("calcSpectrum: Can not calculate %s\n",lop);
            if(errcode) *errcode=-1;
            return 0;
          }
                                            
  for(n=0;n<Nodd;n++) if(strcmp(lop,OddPrtcls[n].name)==0 ||
                         strcmp(lop,OddPrtcls[n].aname)==0  ) break;
  name=OddPrtcls[n].name;
  aname=OddPrtcls[n].aname;
  vcs=calcSpectrum0(name,aname, key, SnmSun,SNmSun,SnmEarth,SNmEarth);

  if(strcmp(name,aname))
  {   
    vcs+=calcSpectrum0(name,name,0,buffnSun,buffNSun,buffnEarth,buffNEarth);
    for(i=0;i<NZ;i++){ SnmSun[i]+= buffnSun[i];  SNmSun[i]+=buffNSun[i];   SnmEarth[i]+= buffnEarth[i];  SNmEarth[i]+=buffNEarth[i];}
  }      

    
  if(vcs) for(i=0;i<NZ;i++){ SnmSun[i]/=vcs; SNmSun[i]/=vcs;SnmEarth[i]/=vcs; SNmEarth[i]/=vcs; }
  if(strcmp(name,aname)) vcs/=2;
  return vcs*2.9979E-26; 
}

static double *SpN_stat=NULL;
static double Enu_stat,Emu_stat,alpha_stat,beta_stat;


static int in2=1; /* proton;   -1 for neutron */
static int in1=1; /* neutrino; -1 for anti-neutrino */


static double dSigmadE_nu2mu(double Enu,double Emu)
{  double GF=1.16637E-5;  /*GeV−2*/
   double mp=0.939;
   double GeVcm=0.50677E14;
   double C=2*GF*GF*mp/M_PI/GeVcm/GeVcm;
   double a=0,b=0;
   switch (in1)
   { case  1:switch (in2)
             { 
               case  1: a=0.15,b=0.04; break;
               case -1: a=0.25,b=0.06; break; 
             }
             break;  
     case -1:switch (in2)
             { 
               case  1: a=0.06,b=0.25; break;
               case -1: a=0.04,b=0.15; break; 
             }
             
             break;
   }
   return C*(a+b*Emu*Emu/(Enu*Enu));       
} 

#define OLD 
#ifdef OLD
static double integrandX(double x)
{   return  dSigmadE_nu2mu(Enu_stat, Emu_stat+x*alpha_stat); }     
#else
static double integrandX(double x)
{  double g=alpha_stat/beta_stat;  
   double ex=exp(x*beta_stat);
   double q= (0.10566 /* muon mass*/)/(65865.4 /* muon life time * c */)/alpha_stat;
   double Emu_prim=(Emu_stat+g)*ex-g;
   return  dSigmadE_nu2mu(Enu_stat, Emu_prim )*ex*pow(Emu_stat*(Emu_stat+g)/(Emu_prim*(Emu_prim+g)),q); 
} 
#endif


#ifdef OLD
static double integrandEnu(double E)
{ 
  Enu_stat=E;
  return  SpectdNdE(Enu_stat,SpN_stat)*simpson(integrandX,0,(Enu_stat-Emu_stat)/alpha_stat,1.E-4);
}
#else
static double integrandEnu(double E)
{ 
  Enu_stat=E;
  double g=alpha_stat/beta_stat;
  double xMax=log((Enu_stat+g)/(Emu_stat+g))/beta_stat;
  return  SpectdNdE(Enu_stat,SpN_stat)*simpson(integrandX,0,xMax,1.E-4);
}
#endif



void nm2mUpward(double * Sp_nm, int inP,  double rho, double P_rate,  double *Sp_mu)
{
   int i;
   double C;
   double  NA=6.002141E23 /* mol-1*/;
 
   Sp_mu[0]=0; 
   
   SpN_stat=Sp_nm;
   in1=inP;  
    
   alpha_stat=0.002*rho;
   beta_stat=3.0E-6*rho;
   for(i=1;i<NZ;i++)
   {  Emu_stat=Mcdm*exp(Zi(i));
      if(Emu_stat<1)  Sp_mu[i]=0; else
      {
         in2=1;
         Sp_mu[i]=simpson(integrandEnu,Emu_stat,Mcdm,1.E-4)*P_rate;
         in2=-1;
         Sp_mu[i]+=simpson(integrandEnu,Emu_stat,Mcdm,1.E-4)*(1-P_rate);
         Sp_mu[i]*=rho*NA*Emu_stat;
      }           
   }                  
}

static double integrandEnuContained(double E) { return  SpectdNdE(E,SpN_stat)*dSigmadE_nu2mu(E,Emu_stat) ; }

void nm2mContained(double * Sp_nm, int inP,  double rho, double P_rate,  double *Sp_mu)
{
   int i;
   double C;
   double  NA=6.002141E23 /* mol-1*/;
 
   Sp_mu[0]=0; 
   
   SpN_stat=Sp_nm;
   in1=inP;  

   for(i=1;i<NZ;i++)
   {  Emu_stat=Mcdm*exp(Zi(i));
 printf("Emu_stat=%E\n",Emu_stat); 
      if(Emu_stat<1)  Sp_mu[i]=0; else
      {
         in2=1;
         Sp_mu[i]=simpson(integrandEnuContained,Emu_stat,Mcdm,1.E-4)*P_rate;
         in2=-1;
         Sp_mu[i]+=simpson(integrandEnuContained,Emu_stat,Mcdm,1.E-4)*(1-P_rate);
         Sp_mu[i]*=rho*NA*Emu_stat*1E5;
      }           
   }                  
}




/*

int main(int argc, char** argv)
{
    int pdg,err;
    double SpNu[NZ],vcs;
    sscanf(argv[1],"%lf",&Mcdm);
    sscanf(argv[2],"%d",&pdg); 
    
    err=basicNuSpectra(pdg,1, SpNu);  
    displaySpectrum(SpNu,"nu spectrum", Mcdm/100,Mcdm,0); 

    err=basicNuSpectra(pdg,-1, SpNu);  
    displaySpectrum(SpNu,"nu-bar  spectrum", Mcdm/100,Mcdm,0); 
  
    err=basicNuSpectra(pdg,0, SpNu);  
    displaySpectrum(SpNu,"(2*nu+1*nu-bar)/3  spectrum",3,Mcdm,0); 
    
    vcs=calcNuSpectrum(4,SpNu,&err);
    displaySpectrum(SpNu," ",3,Mcdm,0);    
}

*/