/*=========   NMSSM scenario  ==========
  One can define SUGRA  for GUT scale scenario
  or EWSB  for low scale input.
  Otherwise the program should read SLHA file.
=======================================*/ 

//#define SUGRA 
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
      
//#define HIGGSBOUNDS "../Packages/HiggsBounds-4.1.0"
//#define HIGGSSIGNALS "../Packages/HiggsSignals-1.2.0"      

#define OMEGA            
      /* Calculate relic density and display contribution of
         individual channels 
      */
//#define INDIRECT_DETECTION  
      /* Compute spectra of gamma/positron/neutrinos
         for DM annihilation; calculate <sigma*v> and
         integrate gamma signal over DM galactic squared
         density for given line of sight.  
      */
//#define LoopGAMMA
      
/*#define RESET_FORMFACTORS  */
      /* Modify default nucleus form factors, 
         DM velocity distribution,
         A-dependence of Fermi-dencity
      */     
#define CDM_NUCLEON     
      /* Calculate amplitudes and cross-sections for 
         CDM-mucleon collisions 
      */  
/*#define CDM_NUCLEUS */ 
      /* Calculate number of events for 1kg*day 
         and recoil energy distibution for various nuclei
      */
#define NEUTRINO 
#define DECAYS
    /* Calculate decay widths and branchings  */

/*#define CROSS_SECTIONS */
      /* Calculate cross sections and widths for 
         reactions specified by the user
      */
/*===== end of Modules  ======*/

/*===== Options ========*/
/* #define SHOWPLOTS */
     /* Display  graphical plots on the screen */ 

#define CLEAN to clean intermediate files


/*===== End of DEFINE  settings ===== */


#include"../sources/micromegas.h"
#include"../sources/micromegas_aux.h"
#include"lib/pmodel.h"


int main(int argc,char** argv)
{  int err,nw;
   char cdmName[10];
   int spin2, charge3,cdim;


  ForceUG=1;   // to Force Unitary Gauge assign 1 
  
#ifdef SUGRA
{
  double m0,mhf,a0,tb;
  double Lambda, aLambda,aKappa,sgn;

  if(argc<7) 
  { 
    printf(" This program needs 6 parameters:\n"
           "   m0      common scalar mass at GUT scale\n"
           "   mhf     common gaugino mass at GUT scale\n"
           "   a0      trilinear soft breaking parameter at GUT scale\n"
           "   tb      tan(beta) \n"
           "   Lambda   Lambda parameter at SUSY\n"
           "   aKappa  aKappa parameter at GUT\n"
           );
    printf(" Auxiliary parameters are:\n"
           "   sgn     +/-1,  sign of Higgsino mass term (default 1)\n" 
           "   aLambda at GUT (default aLambda=a0)\n"    
           "   Mtp     top quark pole mass\n"
           "   MbMb    Mb(Mb) scale independent b-quark mass\n"
           "   alfSMZ  strong coupling at MZ\n");
    printf("Example:  ./main 135 600 -1300 2 0.5 -1400\n");
      exit(1); 
  } else  
  {  double Mtp,MbMb,alfSMZ;
     sscanf(argv[1],"%lf",&m0);
     sscanf(argv[2],"%lf",&mhf);
     sscanf(argv[3],"%lf",&a0);
     sscanf(argv[4],"%lf",&tb);
     sscanf(argv[5],"%lf",&Lambda);
     sscanf(argv[6],"%lf",&aKappa);
     if(argc>7)  sscanf(argv[7],"%lf",&sgn); else sgn=1;
     if(argc>8)  sscanf(argv[8],"%lf",&aLambda); else aLambda=a0;

     if(argc>9){ sscanf(argv[9],"%lf",&Mtp);    assignValW("Mtp",Mtp);      }
     if(argc>10){ sscanf(argv[10],"%lf",&MbMb);   assignValW("MbMb",MbMb);    }
     if(argc>11){ sscanf(argv[11],"%lf",&alfSMZ); assignValW("alfSMZ",alfSMZ);}
  }

  err=nmssmSUGRA( m0,mhf, a0,tb, sgn,  Lambda, aLambda, aKappa); 
}
#elif defined(EWSB)
{
  if(argc!=2)
  { 
      printf(" Correct usage:  ./main <file with NMSSM parameters> \n");
      printf(" Example      :  ./main  data1.par \n");
      exit(1);
  }

  err=readVarNMSSM(argv[1]);
  
  if(err==-1)     {printf("Can not open the file\n"); exit(1);}
  else if(err>0)  { printf("Wrong file contents at line %d\n",err);exit(1);}          
  err=nmssmEWSB(0);  
}
#else
{
   printf("\n========= SLHA file input =========\n");

   if(argc <2) 
   {  printf("The program needs one argument:the name of SLHA input file.\n"
            "Example: ./main spectr.dat \n");
      exit(1);
   }  
   
   printf("Initial file  \"%s\"\n",argv[1]); 
   err=lesHinput(argv[1]);
}

#endif

  slhaWarnings(stdout);
  
  if(err) exit(1);

  err=sortOddParticles(cdmName);


  if(err) { printf("Can't calculate %s\n",cdmName); return 1;}

  err=calcMainFunc();   

  qNumbers(cdmName,&spin2, &charge3, &cdim);
  printf("\nDark matter candidate is '%s' with spin=%d/2\n",
  cdmName,       spin2); 
  if(charge3) { printf("Dark Matter has electric charge %d/3\n",charge3); exit(1);}
  if(cdim!=1) { printf("Dark Matter is a color particle\n"); exit(1);}
  if(strcmp(cdmName,"~o1")) printf(" ~o1 is not CDM\n"); 
                    else o1Contents(stdout);

// If you would like to use widths calculated by NMSSMTools, 
// uncomment the line below

//   slhaRead("decay",1+8);

#ifdef MASSES_INFO
{
  printf("\n=== MASSES OF HIGGS AND SUSY PARTICLES: ===\n");
  printHiggs(stdout);
  printMasses(stdout,1);
}
#endif

#ifdef CONSTRAINTS
{ double constr0,constrM, constrP;

  printf("\n\n==== Physical Constraints: =====\n");

  constr0=bsgnlo(&constrM,&constrP);
  printf("B->s,gamma = %.2E (%.2E ,  %.2E  ) \n",constr0,constrM, constrP );

  constr0= bsmumu(&constrM,&constrP);
  printf("Bs->mu,mu  = %.2E (%.2E ,  %.2E  ) \n",constr0,constrM, constrP );
  
  constr0=btaunu(&constrM,&constrP);
  printf("B+->tau+,nu= %.2E (%.2E ,  %.2E  ) \n",constr0, constrM, constrP );
  
  constr0=deltaMd(&constrM,&constrP);
  printf("deltaMd    = %.2E (%.2E ,  %.2E  ) \n",constr0,constrM, constrP );

  constr0=deltaMs(&constrM,&constrP);
  printf("deltaMs    = %.2E (%.2E ,  %.2E  ) \n",constr0,constrM, constrP );

  constr0=gmuon(&constrM,&constrP);
  printf("(g-2)/BSM = %.2E (%.2E ,  %.2E  ) \n",constr0,constrM, constrP );
     
}
#endif

#ifdef HIGGSBOUNDS
   if(access(HIGGSBOUNDS "/HiggsBounds",X_OK )) system( "cd " HIGGSBOUNDS "; ./configure; make ");
   system("cat spectr decay > HB.in");
   HBblocks("HB.in");
   System("%s/HiggsBounds  LandH SLHA 5 1 HB.in HB.out >hb.stdout",HIGGSBOUNDS);
   slhaRead("HB.out",1+4);
   printf("HB result= %.0E  obsratio=%.2E\n",slhaVal("HiggsBoundsResults",0.,2,1,2), slhaVal("HiggsBoundsResults",0.,2,1,3)  );
   { char hbInfo[100];
    if(0==slhaSTRFormat("HiggsBoundsResults","1 5 ||%[^|]||",hbInfo)) printf("Channel: %s\n",hbInfo);
   }  
#endif

#ifdef HIGGSSIGNALS
#define DataSet " latestresults "
//#define Method  " peak " 
//#define  Method " mass "
#define  Method " both "
#define PDF  " 2 "  // Gaussian
//#define PDF " 1 "  // box 
//#define PDF " 3 "  // box+Gaussia
#define dMh " 2 "
   printf("HiggsSignals:\n");
   if(access(HIGGSSIGNALS "/HiggsSignals",X_OK )) system( "cd " HIGGSSIGNALS "; ./configure; make ");
     system("rm -f HS.in HS.out");
     system("cat spectr decay > HS.in");
     HBblocks("HS.in");
     system("echo 'BLOCK DMASS\n 25 " dMh " '>> HS.in");
     System(HIGGSSIGNALS "/HiggsSignals" DataSet Method  PDF " SLHA 3 1 HS.in > hs.stdout");
     System("grep -A 10000  HiggsSignalsResults HS.in > HS.out");
     slhaRead("HS.out",1+4);
     printf("  Number of observables %.0f\n",slhaVal("HiggsSignalsResults",0.,1,7));
     printf("  total chi^2= %.1E\n",slhaVal("HiggsSignalsResults",0.,1,12));
     printf("  HS p-value = %.1E\n", slhaVal("HiggsSignalsResults",0.,1,13));
#undef dMh
#undef PDF
#undef Method
#undef DataSet
     
#endif


#ifdef OMEGA
{ int fast=1;
  double Beps=1.E-3, cut=0.01;
  double Omega,Xf;
  
// to exclude processes with virtual W/Z in DM   annihilation      
    VZdecay=0; VWdecay=0; cleanDecayTable(); 

// to include processes with virtual W/Z  also  in co-annihilation 
//   VZdecay=2; VWdecay=2; cleanDecayTable(); 
     
  printf("\n==== Calculation of relic density =====\n");  
  Omega=darkOmega(&Xf,fast,Beps);
  printf("Xf=%.2e Omega=%.2e\n",Xf,Omega);
  printChannels(Xf,cut,Beps,1,stdout);
 
  VZdecay=1; VWdecay=1; cleanDecayTable();  // restore default  
}
#endif

#ifdef INDIRECT_DETECTION
{ 
  int err,i;
  double Emin=0.1,/* Energy cut  in GeV   */  sigmaV;
  double vcs_gz,vcs_gg;
  char txt[100];
  double SpA[NZ],SpE[NZ],SpP[NZ];
  double FluxA[NZ],FluxE[NZ],FluxP[NZ];
//  double * SpNe=NULL,*SpNm=NULL,*SpNl=NULL;
  double SpNe[NZ],SpNm[NZ],SpNl[NZ];
  double Etest=Mcdm/2;
  
printf("\n==== Indirect detection =======\n");  

  sigmaV=calcSpectrum(2+4,SpA,SpE,SpP,SpNe,SpNm,SpNl ,&err);
    /* Returns sigma*v in cm^3/sec.     SpX - calculated spectra of annihilation.
       Use SpectdNdE(E, SpX) to calculate energy distribution in  1/GeV units.
       
       First parameter 1-includes W/Z polarization
                       2-includes gammas for 2->2+gamma
                       4-print cross sections             
    */

  if(SpA)
  { 
     double fi=0.1,dfi=0.05; /* angle of sight and 1/2 of cone angle in [rad] */ 

     gammaFluxTab(fi,dfi, sigmaV, SpA,  FluxA);
     printf("\nPhoton flux  for angle of sight f=%.2f[rad]\n"
     "and spherical region described by cone with angle %.2f[rad]\n",fi,2*dfi);
     
#ifdef SHOWPLOTS
     sprintf(txt,"Photon flux[cm^2 s GeV]^{1} at f=%.2f[rad], cone angle %.2f[rad]",fi,2*dfi);
     displaySpectrum(FluxA,txt,Emin,Mcdm,1);
#endif
     printf("Photon flux = %.2E[cm^2 s GeV]^{-1} for E=%.1f[GeV]\n",SpectdNdE(Etest, SpA), Etest);       
  }

  if(SpE)
  { 
    posiFluxTab(Emin, sigmaV, SpE,  FluxE);
#ifdef SHOWPLOTS     
    displaySpectrum(FluxE,"positron flux [cm^2 s sr GeV]^{-1}" ,Emin,Mcdm,1);
#endif
    printf("\nPositron flux  =  %.2E[cm^2 sr s GeV]^{-1} for E=%.1f[GeV] \n",
    SpectdNdE(Etest, FluxE),  Etest);           
  }
  
  if(SpP)
  { 
    pbarFluxTab(Emin, sigmaV, SpP, FluxP  ); 
#ifdef SHOWPLOTS    
     displaySpectrum(FluxP,"antiproton flux [cm^2 s sr GeV]^{-1}" ,Emin,Mcdm,1);
#endif
    printf("\nAntiproton flux  =  %.2E[cm^2 sr s GeV]^{-1} for E=%.1f[GeV] \n",
    SpectdNdE(Etest, FluxP),  Etest);             
  }
}  
#endif

#ifdef LoopGAMMA  
  { double vcs_gg,vcs_gz;
    if(loopGamma(&vcs_gg, &vcs_gz )==0)
    {
      printf("Gamma  ray lines:\n");
      printf("E=%.2E[GeV]  vcs(Z,A)= %.2E[cm^3/s]\n",Mcdm-91.19*91.19/4/Mcdm,vcs_gz);  
      printf("E=%.2E[GeV]  vcs(A,A)= %.2E[cm^3/s]\n",Mcdm,vcs_gg);
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

   calcScalarQuarkFF( Mu/Md, Ms/Md, sigmaPiN[MeV], sigma0[MeV])  
   calculates and rewrites Scalar form factors
*/

  printf("protonFF (default) d %E, u %E, s %E\n",ScalarFFPd, ScalarFFPu,ScalarFFPs);                               
  printf("neutronFF(default) d %E, u %E, s %E\n",ScalarFFNd, ScalarFFNu,ScalarFFNs);

    calcScalarQuarkFF(0.553,18.9,45.5,26.);

//  To restore default form factors of  version 2  call 
//  calcScalarQuarkFF(0.553,18.9,55.,243.5);
 

  printf("protonFF (new)     d %E, u %E, s %E\n",ScalarFFPd, ScalarFFPu,ScalarFFPs);                               
  printf("neutronFF(new)     d %E, u %E, s %E\n",ScalarFFNd, ScalarFFNu,ScalarFFNs);

}
#endif

#ifdef CDM_NUCLEON
{ double pA0[2],pA5[2],nA0[2],nA5[2];
  double Nmass=0.939; /*nucleon mass*/
  double SCcoeff;        

printf("\n==== Calculation of CDM-nucleons amplitudes  =====\n");   

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
 
#ifdef NEUTRINO
{ double nu[NZ], nu_bar[NZ],mu[NZ];
  double Ntot;
  int forSun=1;
  double Emin=0.01;

  printf("\n===============Neutrino Telescope=======  for  ");
  if(forSun) printf("Sun\n"); else printf("Earth\n");

  err=neutrinoFlux(Maxwell,forSun, nu,nu_bar);
#ifdef SHOWPLOTS
  displaySpectrum(nu,"nu flux from Sun [1/Year/km^2/GeV]",Emin,Mcdm,1);
  displaySpectrum(nu_bar,"nu-bar from Sun [1/Year/km^2/GeV]",Emin,Mcdm,1);
#endif
{ double Ntot;
  double Emin=10; //GeV
  spectrInfo(Emin/Mcdm,nu, &Ntot,NULL);
    printf(" E>%.1E GeV neutrino flux       %.3E [1/Year/km^2] \n",Emin,Ntot);
  spectrInfo(Emin/Mcdm,nu_bar, &Ntot,NULL);
    printf(" E>%.1E GeV anti-neutrino flux  %.3E [1/Year/km^2]\n",Emin,Ntot);
}

/* Upward events */

  muonUpward(nu,nu_bar, mu);
#ifdef SHOWPLOTS  
  displaySpectrum(mu,"Upward muons[1/Year/km^2/GeV]",1,Mcdm/2,1);
#endif
  { double Ntot;
    double Emin=1; //GeV
    spectrInfo(Emin/Mcdm,mu, &Ntot,NULL);
    printf(" E>%.1E GeV Upward muon flux    %.3E [1/Year/km^2]\n",Emin,Ntot);
  }

/* Contained events */
  muonContained(nu,nu_bar,1., mu);
#ifdef SHOWPLOTS  
  displaySpectrum(mu,"Contained  muons[1/Year/km^3/GeV]",Emin,Mcdm,1);
#endif
  { double Ntot;
    double Emin=1; //GeV
    spectrInfo(Emin/Mcdm,mu, &Ntot,NULL);
    printf(" E>%.1E GeV Contained muon flux %.3E [1/Year/km^3]\n",Emin,Ntot);
  }
}
#endif

#ifdef DECAYS
{  
  txtList L;
   double width,br;
   char * pname;
   
   printf("\nParticle decays\n"); 
   pname = "h1";
    width=pWidth(pname,&L);
    printf("%s :   total width=%.3E \n and Branchings:\n",pname,width);
    printTxtList(L,stdout);

   printf("\nParticle decays\n");
   pname = "ha";
   width=pWidth(pname,&L);
   printf("%s-> :   total width=%.3E \n and Branchings:\n",pname,width);
   printTxtList(L,stdout);

   pname = "~o2";
    width=pWidth(pname,&L);
    printf("%s:   total width=%.3E \n and Branchings:\n",pname,width);
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
  cc=newProcess("e%,E%->2*x");
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

#ifdef CLEAN 
  system("rm -f inp spectr decay omega nngg.out HB.slha  Key.dat");
  system("rm -f HB.in HB.out HS.in HS.out hb.stdout hs.stdout debug_channels.txt debug_predratio.txt");   
#endif 

  return 0;
}

