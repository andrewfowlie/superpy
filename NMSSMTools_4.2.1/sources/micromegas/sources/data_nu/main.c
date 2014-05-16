/*======  Spectrum calculator  ========= 
   Choose RGE from the list below. SuSpect is included 
   in micrOMEGAs, to use another code define the path 
   to the corresponding package in lib/Makefile
=====================================*/ 
#define RGE suspect
     /* choose 'suspect','isajet','softSusy','spheno'*/

/*=========   SUSY scenario  ==========
  One can define SUGRA, AMSB, EWSB (for low scale input). 
  By default program reads SLHA data file 
=======================================*/
//#define SUGRA 
// #define AMSB  
#define EWSB 

/*====== Modules ===============
   Keys to switch on 
   various modules of micrOMEGAs
================================*/

#define MASSES_INFO      
      /* Display information about SUSY and Higgs masses 
      */
#define CONSTRAINTS     
      /* Display  deltarho, B_>sgamma, Bs->mumu, gmuon and
         check LEP mass limits 
      */ 
#define OMEGA            
      /* Calculate relic density and display contribution of
         individual channels 
      */
#define INDIRECT_DETECTION  
      /* Compute spectra of gamma/positron/neutrinos
         for DM annihilation; calculate <sigma*v> and
         integrate gamma signal over DM galactic squared
         density for given line of sight.  
      */
      
/*#define RESET_FORMFACTORS*/
      /* Modify default nucleus form factors, 
         DM velocity distribution,
         A-dependence of Fermi-dencity
      */     
#define CDM_NUCLEON 
      /* Calculate amplitudes and cross-sections for 
         CDM-mucleon collisions 
      */  
/* #define TEST_Direct_Detection */
      /* 
        Compare analytical formula for DD against micrOMEGAS calculation.
        As well as compare tree level and box improved approaches.
       */      
/*#define CDM_NUCLEUS*/ 
      /* Calculate number of events for 1kg*day 
         and recoil energy distibution for various nuclei
      */
/*#define DECAYS */
      /* Calculate decay widths and branchings  */      
/* #define CROSS_SECTIONS */
      /* Calculate cross sections of reactions specified by the user */

/*===== end of Modules  ======*/

/*===== Options ========*/
//#define SHOWPLOTS 
     /* Display  graphical plots on the screen */ 

/*===== End of DEFINE  settings ===== */


#include"../sources/micromegas.h"
#include"../sources/micromegas_aux.h"
#include"lib/pmodel.h"

#include"nuRead.c"
#include"sun.c"
#include"earth.c"

#define SUGRAMODEL_(A) A ## SUGRA
#define SUGRAMODEL(A) SUGRAMODEL_(A)

#define AMSBMODEL_(A) A ## AMSB
#define AMSBMODEL(A) AMSBMODEL_(A)

#define EWSBMODEL_(A) A ## EwsbMSSM
#define EWSBMODEL(A) EWSBMODEL_(A)

#define PRINTRGE_(A) printf(" Spectrum calculator is %s\n", #A)
#define PRINTRGE(A)  PRINTRGE_(A)


int main(int argc,char** argv)
{  int err;
   char cdmName[10];
   int spin2, charge3,cdim;
   double * massPtr;

 delFiles=1; /* switch to save/delete RGE input/output */
 ForceUG=0;  /* to Force Unitary Gauge assign 1 */
#ifdef SUGRA
{
  double m0,mhf,a0,tb;
  double gMG1, gMG2, gMG3,  gAl, gAt, gAb,  sgn, gMHu,  gMHd,
         gMl2, gMl3, gMr2, gMr3, gMq2, gMq3, gMu2, gMu3, gMd2, gMd3;
         
  printf("\n========= mSUGRA scenario =====\n");
  PRINTRGE(RGE);

  if(argc<5) 
  { 
    printf(" This program needs 4 parameters:\n"
           "   m0      common scalar mass at GUT scale\n"
           "   mhf     common gaugino mass at GUT scale\n"
           "   a0      trilinear soft breaking parameter at GUT scale\n"
           "   tb      tan(beta) \n");
    printf(" Auxiliary parameters are:\n"
           "   sgn     +/-1,  sign of Higgsino mass term (default 1)\n"    
           "   Mtp     top quark pole mass\n"
           "   MbMb    Mb(Mb) scale independent b-quark mass\n"
           "   alfSMZ  strong coupling at MZ\n");
/*    printf("Example: ./main 70 250 -300 10\n");  */
      printf("Example: ./main 120 500 -350 10 1 173.1 \n");
      exit(1); 
  } else  
  {  double Mtp,MbMb,alfSMZ;
     sscanf(argv[1],"%lf",&m0);
     sscanf(argv[2],"%lf",&mhf);
     sscanf(argv[3],"%lf",&a0);
     sscanf(argv[4],"%lf",&tb);
     if(argc>5)sscanf(argv[5],"%lf",&sgn); else sgn=1;
     if(argc>6){ sscanf(argv[6],"%lf",&Mtp);    assignValW("Mtp",Mtp);      }
     if(argc>7){ sscanf(argv[7],"%lf",&MbMb);   assignValW("MbMb",MbMb);    }
     if(argc>8){ sscanf(argv[8],"%lf",&alfSMZ); assignValW("alfSMZ",alfSMZ);}
  }

/*==== simulation of mSUGRA =====*/
  gMG1=mhf, gMG2=mhf,gMG3=mhf;
  gAl=a0,   gAt=a0,  gAb=a0;  gMHu=m0,  gMHd=m0;
  gMl2=m0,  gMl3=m0, gMr2=m0, gMr3=m0;
  gMq2=m0,  gMq3=m0, gMu2=m0, gMd2=m0, gMu3=m0, gMd3=m0;

  err= SUGRAMODEL(RGE) (tb,  
    gMG1, gMG2, gMG3,  gAl,  gAt, gAb,  sgn, gMHu, gMHd,
    gMl2, gMl3, gMr2, gMr3, gMq2,  gMq3, gMu2, gMu3, gMd2, gMd3); 
}
#elif defined(AMSB)
{
  double m0,m32,sgn,tb;

  printf("\n========= AMSB scenario =====\n");
  PRINTRGE(RGE);
  if(argc<4) 
  { 
    printf(" This program needs 3 parameters:\n"
           "   m0      common scalar mass at GUT scale\n"
           "   m3/2    gravitino mass\n"
           "   tb      tan(beta) \n");
    printf(" Auxiliary parameters are:\n"
           "   sgn     +/-1,  sign of Higgsino mass term (default 1)\n"    
           "   Mtp     top quark pole mass\n"
           "   MbMb    Mb(Mb) scale independent b-quark mass\n"
           "   alfSMZ  strong coupling at MZ\n");
   printf("Example: ./main 450  60000 10\n");                                                                          
   exit(1); 
  } else  
  {  double Mtp,MbMb,alfSMZ;
     sscanf(argv[1],"%lf",&m0);
     sscanf(argv[2],"%lf",&m32);
     sscanf(argv[3],"%lf",&tb);
     if(argc>4)sscanf(argv[4],"%lf",&sgn); else sgn=1;
     if(argc>5){ sscanf(argv[5],"%lf",&Mtp);    assignValW("Mtp",Mtp);      }
     if(argc>6){ sscanf(argv[6],"%lf",&MbMb);   assignValW("MbMb",MbMb);    }
     if(argc>7){ sscanf(argv[7],"%lf",&alfSMZ); assignValW("alfSMZ",alfSMZ);}
  }

  err= AMSBMODEL(RGE)(m0,m32,tb,sgn);
 
}
#elif defined(EWSB)
{ 
   printf("\n========= EWSB scale input =========\n");
   PRINTRGE(RGE);

   if(argc <2) 
   {  printf("The program needs one argument:the name of file with MSSM parameters.\n"
            "Example: ./main mssm1.par \n");
      exit(1);
   }  
   
   printf("Initial file  \"%s\"\n",argv[1]);
     
   err=readVarMSSM(argv[1]);
          
   if(err==-1)     { printf("Can not open the file\n"); exit(2);}
   else if(err>0)  { printf("Wrong file contents at line %d\n",err);exit(3);}

   err=EWSBMODEL(RGE)();
}
#else
{
   printf("\n========= SLHA file input =========\n");

   if(argc <2) 
   {  printf("The program needs one argument:the name of SLHA input file.\n"
            "Example: ./main suspect2_lha.out \n");
      exit(1);
   }  
   
   printf("Initial file  \"%s\"\n",argv[1]);
   err=lesHinput(argv[1]);
   if(err) exit(2);
}
#endif
  
  { int nw;
    printf("Warnings from spectrum calculator:\n");
    nw=slhaWarnings(stdout);
    if(nw==0) printf(" .....none\n");
  }  
  if(err) exit(1);
  err=sortOddParticles(cdmName);
  if(err) { printf("Can't calculate %s\n",cdmName); return 1;}

  qNumbers(cdmName,&spin2, &charge3, &cdim);
  printf("\nDark matter candidate is '%s' with spin=%d/2  mass=%.2E\n",
  cdmName,       spin2, Mcdm); 
    
  
  if(charge3) { printf("Dark Matter has electric charge %d/3\n",charge3); exit(1);}
  if(cdim!=1) { printf("Dark Matter is a color particle\n"); exit(1);}
  if(strcmp(cdmName,"~o1")) printf(" ~o1 is not CDM\n"); 
                              else o1Contents(stdout);

                
#ifdef MASSES_INFO
{
  printf("\n=== MASSES OF HIGGS AND SUSY PARTICLES: ===\n");
  printHiggs(stdout);
  printMasses(stdout,1);
}
#endif

#ifdef CONSTRAINTS
{ printf("\n\n==== Physical Constraints: =====\n"); 
  printf("deltartho=%.2E\n",deltarho());
  printf("gmuon=%.2E\n", gmuon());
  printf("bsgnlo=%.2E\n", bsgnlo());
  printf("bsmumu=%.2E\n", bsmumu());
  printf("btaunu=%.2E\n", btaunu());
  if(masslimits()==0) printf("MassLimits OK\n");
}
#endif

#ifdef OMEGA
{ int fast=1;
  double Beps=1.E-5, cut=0.01;
  double Omega,Xf;   
  printf("\n==== Calculation of relic density =====\n");  
  Omega=darkOmega(&Xf,fast,Beps);
  printf("Xf=%.2e Omega=%.2e\n",Xf,Omega);
  printChannels(Xf,cut,Beps,1,stdout);
}
#endif


#ifdef INDIRECT_DETECTION
{ 
  int err,i;
  double Emin=1,SMmev=320;/*Energy cut in GeV and solar potential in MV*/
  double  sigmaV;
  double vcs_gz,vcs_gg;
  char txt[100];
  double SpA[NZ],SpE[NZ],SpP[NZ];
  double * SpNe=NULL,*SpNm=NULL,*SpNl=NULL;
  double Etest=Mcdm/2;
 
/* default DarkSUSY parameters */

/*
    K_dif=0.036;
    L_dif=4;  
    Delta_dif=0.6; 
    Vc_dif=10;
    Rdisk=30;
    SMmev=320;
*/                        
  
printf("\n==== Indirect detection =======\n");  

  sigmaV=calcSpectrum( 1+2+4,SpA,SpE,SpP,SpNe,SpNm,SpNl ,&err);
    /* Returns sigma*v in cm^3/sec.     SpX - calculated spectra of annihilation.
       Use SpectdNdE(E, SpX) to calculate energy distribution in  1/GeV units.
       
       First parameter 1-includes W/Z polarization
                       2-includes gammas for 2->2+gamma
                       4-print cross sections             
    */
  printf("sigmav=%.2E[cm^3/s]\n",sigmaV);  

  if(SpA)
  { 
     double fi=0.,dfi=M_PI/180.; /* angle of sight and 1/2 of cone angle in [rad] */ 
                                                   /* dfi corresponds to solid angle 1.E-3sr */                                             
     printf("Photon flux  for angle of sight f=%.2f[rad]\n"
     "and spherical region described by cone with angle %.4f[rad]\n",fi,2*M_PI*(1-cos(dfi)));
     gammaFluxTab(fi,dfi, sigmaV, SpA, SpA);

#ifdef SHOWPLOTS
     sprintf(txt,"Photon flux for angle of sight %.2f[rad] and cone angle %.2f[rad]",fi,2*dfi);
     displaySpectrum(SpA,txt,Emin,0.98*Mcdm,1);
#endif
     printf("Photon flux = %.2E[cm^2 s GeV]^{-1} for E=%.1f[GeV]\n",SpectdNdE(Etest, SpA), Etest);       
     if(loopGamma(&vcs_gz,&vcs_gg)==0)
     {
         printf("Gamma  ray lines:\n");
         printf("E=%.2E[GeV]  vcs(Z,A)= %.2E[cm^3/s], flux=%.2E[cm^2 s]^{-1}\n",Mcdm-91.19*91.19/4/Mcdm,vcs_gz,
                               gammaFlux(fi,dfi,vcs_gz));  
         printf("E=%.2E[GeV]  vcs(A,A)= %.2E[cm^3/s], flux=%.2E[cm^2 s]^{-1}\n",Mcdm,vcs_gg, 
                             2*gammaFlux(fi,dfi,vcs_gg));
     }
  }

  if(SpE)
  { 

    posiFluxTab(Emin, sigmaV, SpE, SpE);
    if(SMmev>0)  solarModulation(SMmev,0.0005,SpE,SpE);    
#ifdef SHOWPLOTS     
    displaySpectrum(SpE,"positron flux [cm^2 s sr GeV]^{-1}" ,Emin,0.95*Mcdm,1);
#endif
    printf("Positron flux  =  %.2E[cm^2 sr s GeV]^{-1} for E=%.1f[GeV] \n",
    SpectdNdE(Etest, SpE),  Etest); 
  }
  
  if(SpP)
  {
    pbarFluxTab(Emin, sigmaV, SpP,  SpP); 
    if(SMmev>0)  solarModulation(SMmev,1,SpP,SpP);     
#ifdef SHOWPLOTS    
     displaySpectrum(SpP,"antiproton flux [cm^2 s sr GeV]^{-1}" ,Emin,0.5*Mcdm,1);
#endif
    printf("Antiproton flux  =  %.2E[cm^2 sr s GeV]^{-1} for E=%.1f[GeV] \n",
    SpectdNdE(Etest, SpP),  Etest);     
  }
}  
#endif

#ifdef RESET_FORMFACTORS
{
/* 
   The user has approach to form factors  which specifies quark contents 
   of  proton and nucleon via global parametes like
      <Type>FF<Nucleon><q>
   where <Type> can be "Scalar", "pVector", and "Sigma"; 
         <Nucleon>     "P" or "N" for proton and neutron
         <q>            "d", "u","s"

   calcScalarFF( Mu/Md, Ms/Md, sigmaPiN[MeV], sigma0[MeV])  
   calculates and rewrites Scalar form factors
*/

  printf("protonFF (default) d %E, u %E, s %E\n",ScalarFFPd, ScalarFFPu,ScalarFFPs);                               
  printf("neutronFF(default) d %E, u %E, s %E\n",ScalarFFNd, ScalarFFNu,ScalarFFNs);

  calcScalarFF(0.553,18.9,70.,35.);

  printf("protonFF (new)     d %E, u %E, s %E\n",ScalarFFPd, ScalarFFPu,ScalarFFPs);                               
  printf("neutronFF(new)     d %E, u %E, s %E\n",ScalarFFNd, ScalarFFNu,ScalarFFNs);



/* Option to change parameters of DM velocity  distribution  */   
   SetfMaxwell(220.,600.);
/* 
    dN  ~  exp(-v^2/arg1^2)*Theta(v-arg2)  d^3v     
    Earth velocity with respect to Galaxy defined by 'Vearth' parameter.
    All parameters are  in [km/s] units.       
*/
}
#endif

#ifdef CDM_NUCLEON
{ double pA0[2],pA5[2],nA0[2],nA5[2];
  double Nmass=0.939; /*nucleon mass*/
  double SCcoeff;        

printf("\n==== Calculation of CDM-nucleons amplitudes  =====\n");   
#ifdef TEST_Direct_Detection
printf("         TREE LEVEL\n");

    MSSMDDtest(0, pA0,pA5,nA0,nA5);
    printf("Analitic formulae\n");
    printf("proton:  SI  %.3E  SD  %.3E\n",pA0[0],pA5[0]);
    printf("neutron: SI  %.3E  SD  %.3E\n",nA0[0],nA5[0]); 

    nucleonAmplitudes(NULL, pA0,pA5,nA0,nA5);
    printf("CDM-nucleon micrOMEGAs amplitudes:\n");
    printf("proton:  SI  %.3E  SD  %.3E\n",pA0[0],pA5[0]);
    printf("neutron: SI  %.3E  SD  %.3E\n",nA0[0],nA5[0]); 


printf("         BOX DIAGRAMS\n");  

    MSSMDDtest(1, pA0,pA5,nA0,nA5);
    printf("Analitic formulae\n");
    printf("proton:  SI  %.3E  SD  %.3E\n",pA0[0],pA5[0]);
    printf("neutron: SI  %.3E  SD  %.3E\n",nA0[0],nA5[0]); 

    
#endif

    nucleonAmplitudes(FeScLoop, pA0,pA5,nA0,nA5);
    printf("CDM-nucleon micrOMEGAs amplitudes:\n");
    printf("proton:  SI  %.3E  SD  %.3E\n",pA0[0],pA5[0]);
    printf("neutron: SI  %.3E  SD  %.3E\n",nA0[0],nA5[0]); 

  SCcoeff=4/M_PI*3.8937966E8*pow(Nmass*Mcdm/(Nmass+ Mcdm),2.);
    printf("CDM-nucleon cross sections[pb]:\n");
    printf(" proton  SI %.3E  SD %.3E\n",SCcoeff*pA0[0]*pA0[0],3*SCcoeff*pA5[0]*pA5[0]);
    printf(" neutron SI %.3E  SD %.3E\n",SCcoeff*nA0[0]*nA0[0],3*SCcoeff*nA5[0]*nA5[0]);
}
#endif
  
#ifdef CDM_NUCLEUS
{ double dNdE[300];
  double nEvents;

printf("\n======== Direct Detection ========\n");    

  nEvents=nucleusRecoil(Maxwell,73,Z_Ge,J_Ge73,S00Ge73,S01Ge73,S11Ge73,FeScLoop,dNdE);

  printf("73Ge: Total number of events=%.2E /day/kg\n",nEvents);
  printf("Number of events in 10 - 50 KeV region=%.2E /day/kg\n",
                                   cutRecoilResult(dNdE,10,50));
                                                                                                         
#ifdef SHOWPLOTS
    displayRecoilPlot(dNdE,"Distribution of recoil energy of 73Ge",0,199);
#endif

  nEvents=nucleusRecoil(Maxwell,131,Z_Xe,J_Xe131,S00Xe131,S01Xe131,S11Xe131,FeScLoop,dNdE);

  printf("131Xe: Total number of events=%.2E /day/kg\n",nEvents);
  printf("Number of events in 10 - 50 KeV region=%.2E /day/kg\n",
                                   cutRecoilResult(dNdE,10,50));                                   
#ifdef SHOWPLOTS
    displayRecoilPlot(dNdE,"Distribution of recoil energy of 131Xe",0,199);
#endif

  nEvents=nucleusRecoil(Maxwell,23,Z_Na,J_Na23,S00Na23,S01Na23,S11Na23,FeScLoop,dNdE);

  printf("23Na: Total number of events=%.2E /day/kg\n",nEvents);
  printf("Number of events in 10 - 50 KeV region=%.2E /day/kg\n",
                                   cutRecoilResult(dNdE,10,50));                                   
#ifdef SHOWPLOTS
    displayRecoilPlot(dNdE,"Distribution of recoil energy of 23Na",0,199);
#endif

  nEvents=nucleusRecoil(Maxwell,127,Z_I,J_I127,S00I127,S01I127,S11I127,FeScLoop,dNdE);

  printf("I127: Total number of events=%.2E /day/kg\n",nEvents);
  printf("Number of events in 10 - 50 KeV region=%.2E /day/kg\n",
                                   cutRecoilResult(dNdE,10,50));                                   
#ifdef SHOWPLOTS
    displayRecoilPlot(dNdE,"Distribution of recoil energy of 127I",0,199);
#endif
  
}
#endif 

#ifdef DECAYS
{  
  txtList L;
   int dim;
   double width,br;
   char * pname;

   pname = "h";
    width=pWidth(pname,&L,&dim);
    printf("%s->%d*x :   total width=%E \n and Branchings:\n",pname,dim,width);
    printTxtList(L,stdout);

   pname = "l";
    width=pWidth(pname,&L,&dim);
    printf("%s->%d*x :   total width=%E \n and Branchings:\n",pname,dim,width);
    printTxtList(L,stdout);
    printf("Br(e,Ne,nl)= %E\n",findBr(L,"e,Ne,nl"));

   pname = "~o2";
    width=pWidth(pname,&L,&dim);
    printf("%s->%d*x :   total width=%E \n and Branchings:\n",pname,dim,width);
    printTxtList(L,stdout);
    
   pname = "~g";
    width=pWidth(pname,&L,&dim);
    printf("%s->%d*x :   total width=%E \n and Branchings:\n",pname,dim,width);
    printTxtList(L,stdout);
    
    
}
#endif

#ifdef CROSS_SECTIONS
{
  double Pcm=500, cosmin=-0.99, cosmax=0.99, cs;
  numout* cc;
printf("\n====== Calculation of cross section ====\n");  

printf(" e^+, e^- annihilation\n");
  Pcm=500.;
  Helicity[0]=0.5;    /* helicity : spin projection on direction of motion   */    
  Helicity[1]=-0.5;   /* helicities ={ 0.5, -0.5} corresponds to vector state */
  printf("Process e,E->2*x at Pcm=%.3E GeV\n",Pcm);
  cc=newProcess("e%,E%->2*x","eE_2x");
  if(cc)
  { int ntot,l;
    char * name[4];
    procInfo1(cc,&ntot,NULL,NULL);
    for(l=1;l<=ntot; l++)
    { int err;
      double cs;
      char txt[100];
      procInfo2(cc,l,name,NULL);
      sprintf(txt,"%3s,%3s -> %3s %3s  ",name[0],name[1],name[2],name[3]);
      cs= cs22(cc,l,Pcm,cosmin,cosmax,&err);
      if(err) printf("%-20.20s    Error\n",txt);
      else if(cs) printf("%-20.20s  %.2E [pb]\n",txt,cs); 
    }
  } 
}

#endif
{ int pdg=15,i;
  double vcs=0,SpnuSun[NZ],SpNuSun[NZ],SpnuEarth[NZ],SpNuEarth[NZ],
         SpmuSunContained[NZ], SpMuSunContained[NZ],C,CE,SCcoeff;
  double SpmuSun[NZ],SpMuSun[NZ],SpmuEarth[NZ],SpMuEarth[NZ];
  double csSIp,csSIn,csSDp,csSDn;
  double vRotation=270,Vmax=600;
  double AU=150E6; /*km*/
  double Rearth=6378.1; /*km*/
  double pA0[2],nA0[5],pA5[2],nA5[2];
  double Nmass=0.940;
  double Ntot;
/*  
    err=basicNuSpectra(15,0, SpNu);  
    displaySpectrum(SpNu,"nu spectrum of tau", Mcdm/100,Mcdm,0); 

    err=basicNuSpectra(6,0, SpNu);  
    displaySpectrum(SpNu,"nu-bar  spectrum of t", Mcdm/100,Mcdm,0); 
  
    err=basicNuSpectra(5,0, SpNu);  
    displaySpectrum(SpNu,"nu  spectrum of b",3,Mcdm,0); 
*/    
/*
     basicNuSpectraSun(24,1,SpnuSun);
     displaySpectrum(SpnuSun,"W->nu (Sun)",3,Mcdm,0);
     basicNuSpectraEarth(24,1,SpnuEarth);
     displaySpectrum(SpnuEarth,"W->nu (Earth)",3,Mcdm,0);
*/
    vcs=calcNuSpectrum(4,SpnuSun,SpNuSun,SpnuEarth,SpNuEarth,&err);
    spectrInfo(1/Mcdm,SpNuSun,&Ntot,NULL);    
    displaySpectrum(SpNuSun,"Nu-bar (from Sun) ",3,Mcdm,0);
    displaySpectrum(SpNuEarth,"Nu-bar(from Earth) ",3,Mcdm,0);
    
    printf("vcs=%.2E  nuTot=%.2E\n",vcs,Ntot);   
    
    SetfMaxwell(vRotation,Vmax); 
    nucleonAmplitudes(FeScLoop, pA0,pA5,nA0,nA5);
    printf("CDM-nucleon micrOMEGAs amplitudes:\n");
    printf("proton:  SI  %.3E  SD  %.3E\n",pA0[0],pA5[0]);
    printf("neutron: SI  %.3E  SD  %.3E\n",nA0[0],nA5[0]); 

    SCcoeff=4/M_PI*3.8937966E8*pow(Nmass*Mcdm/(Nmass+ Mcdm),2.);
    csSIp=SCcoeff*pA0[0]*pA0[0];
    csSIn=SCcoeff*nA0[0]*nA0[0];
    csSDp=3*SCcoeff*pA5[0]*pA5[0];
    csSDn=3*SCcoeff*nA5[0]*nA5[0];
printf(" csSIp=%.2E  csSIn=%.2E csSDp=%.2E csSDn=%.2E\n",csSIp,csSIn,csSDp,csSDn);
    SetfMaxwell(vRotation,Vmax);    
    C=capture(Mcdm, csSIp,csSIn,csSDp,csSDn,vRotation );
    CE=captureE(Mcdm, csSIp, csSIn, csSDp,csSDn,vRotation);
    printf("Capture %E(Sun)  %E(Earth) \n",C,CE);
        
    for(i=0;i<NZ;i++) {SpNuSun[i]*=C/2/AU/AU/4/M_PI*31556925.2;
                       SpnuSun[i]*=C/2/AU/AU/4/M_PI*31556925.2;
                       SpNuEarth[i]*=CE/2/Rearth/Rearth/4/M_PI*31556925.2;
                       SpnuEarth[i]*=CE/2/Rearth/Rearth/4/M_PI*31556925.2;   

                       }
/* 
    displaySpectrum(SpNuSun,"anti-nu at Earth in /s/km^2",1,Mcdm,0); 
    displaySpectrum(SpnuSun,"    -nu at Earth in /s/km^2",1,Mcdm,0);   
*/
    spectrInfo(1/Mcdm,SpnuSun,&Ntot,NULL);
    printf("Flux_nu(Sun)=%.2E[1/km^2/y]\n",Ntot); 

    spectrInfo(1/Mcdm,SpNuSun,&Ntot,NULL);
    printf("Flaux_NU(Sun)=%.2E[1/km^2/y]\n",Ntot); 
     
    spectrInfo(1/Mcdm,SpnuEarth,&Ntot,NULL);
    printf("Flux_nu(Earth)=%.2E[1/km^2/y]\n",Ntot); 

    spectrInfo(1/Mcdm,SpNuEarth,&Ntot,NULL);
    printf("Flaux_NU(Earth)=%.2E[1/km^2/y]\n",Ntot); 
     
     
     
     
      
    nm2mUpward(SpnuSun,1,2.6, 0.5,  SpmuSun);
    displaySpectrum(SpmuSun," muon Upward flux in /y/km^2",1,Mcdm,1);    
      
    nm2mUpward(SpNuSun,-1,2.6, 0.5,  SpMuSun);
    displaySpectrum(SpMuSun," anti-muon Upward  flux in /y/km^2",1,Mcdm,1);

    nm2mContained(SpnuSun,1,2.6, 0.5,  SpmuSunContained);
    displaySpectrum(SpmuSunContained," muon  Contained  in /y/km^3",1,Mcdm,1);    
      
    nm2mContained(SpNuSun,-1,2.6,0,  SpMuSunContained);
    displaySpectrum(SpMuSunContained," anti-muon Contained  in /y/km^3",1,Mcdm,1);
    
    
    spectrInfo(1/Mcdm,SpmuSun,&Ntot,NULL);
    printf("Flux_mu=%.2E[1/km^2/y\n",Ntot); 

    spectrInfo(1/Mcdm,SpMuSun,&Ntot,NULL);
    printf("Flaux_Mu=%.2E[1/km^2/y\n",Ntot);          
}
  return 0;
}
