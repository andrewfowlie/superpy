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
//#define HIGGSBOUNDS  "../Packages/HiggsBounds-4.1.0"
//#define HIGGSSIGNALS "../Packages/HiggsSignals-1.2.0"

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

//#define LoopGAMMA
      /* Calculate discrete  photon spectrum caused by annihilation of 
         neutralinos into two photons and Z-photon
      */ 

//#define RESET_FORMFACTORS 
      /* Modify default nucleus form factors */     
#define CDM_NUCLEON 
      /* Calculate amplitudes and cross-sections for 
         CDM-mucleon collisions 
      */  
/*#define CDM_NUCLEUS */     
      /* Calculate number of events for 1kg*day 
         and recoil energy distibution for various nuclei
      */
#define NEUTRINO // neutrino telescope signals 
      
#define DECAYS
    /* Calculate decay widths and branchings  */
      
/* #define CROSS_SECTIONS */
      /* Calculate cross sections and widths for 
         reactions specified by the user
      */
/*===== end of Modules  ======*/

/*===== Options ========*/
/*#define SHOWPLOTS */
     /* Display  graphical plots on the screen */ 

#define CLEAN   to clean intermidiate files 

/*===== End of DEFINE  settings ===== */


#include"../sources/micromegas.h"
#include"../sources/micromegas_aux.h"
#include"lib/pmodel.h"


int main(int argc,char** argv)
{  int err;
   char cdmName[10];
   int spin2, charge3,cdim;

 ForceUG=0;    /* to Force Unitary Gauge assign 1 */

  if(argc==1)
  { 
      printf(" Correct usage:  ./main  <file with NMSSM parameters> \n");
      printf("Example: ./main data1.par\n");
      exit(1);
  }
                               
/*  err=readVar(argv[1]);*/
  err=readVarCPVMSSM(argv[1]);
  
  if(err==-1)     {printf("Can not open the file\n"); exit(1);}
  else if(err>0)  { printf("Wrong file contents at line %d\n",err);exit(1);}
           

  if(err) 
  { printf("Problem with RGE solution or spectrum calculation\n");
     exit(1);
  }

  err=sortOddParticles(cdmName);
  if(err) { printf("Can't calculate %s\n",cdmName); return 1;}
   qNumbers(cdmName,&spin2, &charge3, &cdim);
   printf("\nDark matter candidate is '%s' with spin=%d/2 \n",
    cdmName,       spin2); 
   if(charge3) { printf("Dark Matter has electric charge %d/3\n",charge3); exit(1);}
   if(cdim!=1) { printf("Dark Matter is a color particle\n"); exit(1);}
   if(strcmp(cdmName,"~o1")) printf(" ~o1 is not CDM\n"); 
                    else o1Contents(stdout);

// If you would like to use widths calculated by NMSSMTools, 
// uncomment the line below

//   slhaRead("cpsuperh2_slha.out",1+8);


#ifdef MASSES_INFO
{
  printf("\n=== MASSES OF HIGGS AND SUSY PARTICLES: ===\n");
  printHiggs(stdout);
  printMasses(stdout,1);
}
#endif

#ifdef CONSTRAINTS
{ double rat;
  FILE *f;
  printf("\n==== Physical Constraints: =====\n");  
  printf("EDM(Thallium)=%.2E\n",               slhaValFormat("FDIPOLE",0.,"1000812050  1 1 %lf")  ) ;
  printf("EDM(Mercury)=%.2E",                  slhaValFormat("FDIPOLE",0.,"1000801990  1 1 %lf"));
  printf("  # %s\n" ,slhaComment) ;
//  printf("EDM(electron)=%.2E\n",               slhaValFormat("dipoleMOM",0.,"1000812030     1 %lf") ) ;

  printf("EDM(muon)=%.2E\n",                   slhaValFormat("FDIPOLE",0.,"13     1 1 %lf"));

  printf("SUSY contr. to muon Mag.Mom.=%.2E\n",slhaValFormat("FDIPOLE",0.,"13     2 1 %lf"));   

  printf("Br(Bs->mu+,mu-)=%.2E\n",             slhaValFormat("FOBS",0.,"531 1 1  %lf  %*f  2  13 -13"));
  printf("Br(Bd->tau+,tau-)=%.2E\n",           slhaValFormat("FOBS",0.,"511 1 1  %lf  %*f  2 -15 15"));
  printf("Br(B->Xs,gamma)=%.2E\n",             slhaValFormat("FOBS",0.,"5   1 1  %lf  %*f  2   3 22"));
  printf("CP_asymm(B->Xs,gamma)=%.2E[%%]\n",   slhaValFormat("FOBS",0.,"5   3 1  %lf  %*f  2   3 22"));
  printf("Br(Bu->tau,nu)/Br(SM)=%.2E\n",       slhaValFormat("FOBS",0.,"521 2 1  %lf  %*f  2 -15 16"));

//  printf("SUSY mass diff, Bd-bar{Bd} =%.2E[1/ps]\n", slhaValFormat("deltaM",0.,"511 1 %lf"));
//  printf("SUSY mass diff, Bs-bar{Bs} =%.2E[1/ps]\n", slhaValFormat("deltaM",0.,"531 1 %lf"));  
  
#ifdef HIGGSBOUNDS

  if(access(HIGGSBOUNDS "/HiggsBounds",X_OK )) system( "cd " HIGGSBOUNDS "; ./configure; make ");
  System( HIGGSBOUNDS "/HiggsBounds  LandH SLHA 3 1  cpsuperh2_slha.out HB.out");
  slhaRead("HB.out",1+4);
  printf("HB result= %.0E  obsratio=%.2E\n",slhaVal("HiggsBoundsResults",0.,2,1,2), slhaVal("HiggsBoundsResults",0.,2,1,3));
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
     System("cp cpsuperh2_slha.out HS.in"); 
     system("echo 'BLOCK DMASS\n 25 " dMh " '>> HS.in");
     system(HIGGSSIGNALS "/HiggsSignals" DataSet Method  PDF  " SLHA 3 1 HS.in > hs.stdout");
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

}
#endif

#ifdef OMEGA
{ int fast=1;
  double Beps=1.E-4, cut=0.01;
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
  double Emin=1,/* Energy cut  in GeV   */  sigmaV;
  double vcs_gz,vcs_gg;
  char txt[100];
  double SpA[NZ],SpE[NZ],SpP[NZ];
  double FluxA[NZ],FluxE[NZ],FluxP[NZ];
  double * SpNe=NULL,*SpNm=NULL,*SpNl=NULL;
  double Etest=Mcdm/2;
  
printf("\n==== Indirect detection =======\n");  

  sigmaV=calcSpectrum(1+4,SpA,SpE,SpP,SpNe,SpNm,SpNl ,&err);
    /* Returns sigma*v in cm^3/sec.     SpX - calculated spectra of annihilation.
       Use SpectdNdE(E, SpX) to calculate energy distribution in  1/GeV units.
       
       First parameter 1-includes W/Z polarization
                       2-includes gammas for 2->2+gamma
                       4-print cross sections             
    */
  printf("sigmav=%.2E[cm^3/s]\n",sigmaV);  


  if(SpA)
  { 
     double fi=0.1,dfi=0.05; /* angle of sight and 1/2 of cone angle in [rad] */ 

     gammaFluxTab(fi,dfi, sigmaV, SpA, FluxA);     
    printf("Photon flux for angle of sight f=%.2f[rad]\n"
     "and spherical region described by cone with angle %.2f[rad]\n",fi,2*dfi);

#ifdef SHOWPLOTS
     sprintf(txt,"Photon flux[cm^2 s GeV]^{1} at f=%.2f[rad], cone angle %.2f[rad] ",fi,2*dfi);
     displaySpectrum(FluxA,txt,Emin,Mcdm,1);
#endif
     printf("Photon flux = %.2E[cm^2 s GeV]^{-1} for E=%.1f[GeV]\n",SpectdNdE(Etest, FluxA),Etest);       
  }

#ifdef LoopGAMMA
{
     double fi=0.1,dfi=0.05; /* angle of sight and 1/2 of cone angle in [rad] */ 
     double vcs_gz,vcs_gg;
     if(loopGamma(&vcs_gz,&vcs_gg)==0)
     {
         printf("Gamma  ray lines:\n");
         printf("E=%.2E[GeV]  vcs(Z,A)= %.2E[cm^3/s], flux=%.2E[cm^2 s]^{-1}\n",Mcdm-91.19*91.19/4/Mcdm,vcs_gz,
                               gammaFlux(fi,dfi,vcs_gz));  
         printf("E=%.2E[GeV]  vcs(A,A)= %.2E[cm^3/s], flux=%.2E[cm^2 s]^{-1}\n",Mcdm,vcs_gg, 
                             2*gammaFlux(fi,dfi,vcs_gg));
     }
}     
#endif     

  if(SpE)
  { 
    posiFluxTab(Emin, sigmaV, SpE, FluxE);
#ifdef SHOWPLOTS     
    displaySpectrum(FluxE,"positron flux [cm^2 s sr GeV]^{-1}" ,Emin,Mcdm,1);
#endif
    printf("Positron flux  =  %.2E[cm^2 sr s GeV]^{-1} for E=%.1f[GeV] \n",
    SpectdNdE(Etest, FluxE),  Etest);           
  }
  
  if(SpP)
  { 
    pbarFluxTab(Emin, sigmaV, SpP,  FluxP  ); 
#ifdef SHOWPLOTS    
     displaySpectrum(FluxP,"antiproton flux [cm^2 s sr GeV]^{-1}" ,Emin,Mcdm,1);
#endif
    printf("Antiproton flux  =  %.2E[cm^2 sr s GeV]^{-1} for E=%.1f[GeV] \n",
    SpectdNdE(Etest, FluxP),  Etest);             
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

/* setRecoilEnergyGrid(0.7,300); */

  nEvents=nucleusRecoil(Maxwell,73,Z_Ge,J_Ge73,SxxGe73,FeScLoop,dNdE);

  printf("73Ge: Total number of events=%.2E /day/kg\n",nEvents);
  printf("Number of events in 10 - 50 KeV region=%.2E /day/kg\n",
                                   cutRecoilResult(dNdE,10,50));
                                                                                                         
#ifdef SHOWPLOTS
    displayRecoilPlot(dNdE,"Distribution of recoil energy of 73Ge",0,199);
#endif

  nEvents=nucleusRecoil(Maxwell,131,Z_Xe,J_Xe131,SxxXe131,FeScLoop,dNdE);

  printf("131Xe: Total number of events=%.2E /day/kg\n",nEvents);
  printf("Number of events in 10 - 50 KeV region=%.2E /day/kg\n",
                                   cutRecoilResult(dNdE,10,50));                                   
#ifdef SHOWPLOTS
    displayRecoilPlot(dNdE,"Distribution of recoil energy of 131Xe",0,199);
#endif

  nEvents=nucleusRecoil(Maxwell,23,Z_Na,J_Na23,SxxNa23,FeScLoop,dNdE);

  printf("23Na: Total number of events=%.2E /day/kg\n",nEvents);
  printf("Number of events in 10 - 50 KeV region=%.2E /day/kg\n",
                                   cutRecoilResult(dNdE,10,50));                                   
#ifdef SHOWPLOTS
    displayRecoilPlot(dNdE,"Distribution of recoil energy of 23Na",0,199);
#endif

  nEvents=nucleusRecoil(Maxwell,127,Z_I,J_I127,SxxI127,FeScLoop,dNdE);

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
   int dim;
   double width,br;
   char * pname;

   pname = "h1";
    width=pWidth(pname,&L);
    printf("Width(%s)=%.3E \n  Branchings:\n",pname,width);
    printTxtList(L,stdout);

   pname = "l";
    width=pWidth(pname,&L);
    printf("Width(%s)=%.3E \n and Branchings:\n",pname,width);
    printTxtList(L,stdout);
    printf("Br(e,Ne,nl)= %.3E\n",findBr(L,"e,Ne,nl"));

   pname = "~o2";
    width=pWidth(pname,&L);
    printf("Width(%s)=%.3E \n and Branchings:\n",pname,width);
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
 
system(" rm -f CPsuperH.in cpsuperh2_slha.out CPsuperH.out Key.dat nngg.in nngg.out debug_channels.txt debug_predratio.txt"); 
system(" rm -f HB.in HB.out hb.stdout  HS.in HS.out hs.stdout"); 
#endif

  killPlots();
  return 0;
}

