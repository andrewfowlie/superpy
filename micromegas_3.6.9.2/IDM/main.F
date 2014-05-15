
C====== Modules =================
C   Keys to switch on
C   various modules of micrOMEGAs
C================================
   
#define MASSES_INFO
C    Display information about SUSY and Higgs masses
#define CONSTRAINTS
C    Strongly depends on the model. 

C#define HIGGSBOUNDS "../Packages/HiggsBounds-4.1.0"
C#define HIGGSSIGNALS "../Packages/HiggsSignals-1.2.0"


#define OMEGA
C    Calculate relic density and display contribution of
C                            individual channels 
#define INDIRECT_DETECTION
C     Compute spectra of gamma/positron/neutrinos
C     for DM annihilation; calculate <sigma*v> and
C     integrate gamma signal over DM galactic squared
C     density for given line of sight.
C#define RESET_FORMFACTORS
C     Modify default nucleus form factors, 
C     DM velocity distribution,
C     A-dependence of Fermi-dencity.
#define CDM_NUCLEON 
C     Calculate amplitudes and cross-sections for
C     CDM-mucleon collisions 
C#define CDM_NUCLEUS
C     Calculate number of events for 1kg*day
C     and recoil energy distibution for various nuclei

#define NEUTRINO for neutrino telescope

#define DECAYS 
C================ end of Modules  ======

C===== Options ========

C#define SHOWPLOTS
C     Display  graphical plots on the screen
#define CLEAN
C===== End of DEFINE  settings ===== 


      IMPLICIT NONE
C! external functions 
#include "../sources/micromegas.fh"
#include "lib/pmodel.fh"
*      external iargc
*      integer iargc
C! local  Parameters   

      INTEGER fast/1/
      REAL*8  Beps/1.E-4/,cut/0.01/
C Variables      
      INTEGER err
      REAL*8 Omega,OmegaAppr,Xf
      character *100 argv, buff
      character*100 mess
      character*20 cdmName
      integer spin2,charge3,cdim,PDGcdm
      real*8  Emin,sigmaV,fi,dfi,Etest
      real*8 SpA(250),Spe(250),SpP(250) 
      real*8 FluxA(250),FluxE(250),FluxP(250) 
C for CDM_NUCLEON
      real*8  pA0(2),pA5(2),nA0(2),nA5(2),Nmass,SCcoeff
C  for CDM_NUCLEUS
      real*8  dNdE(200),rho,nEvents
C  for NEUTRINO
      integer forSun
      real*8 IntFlux, rhoDet   
      real*8 FluxN(NZ),FluxNb(NZ),FluxMu(NZ)
C  For DECAYS
      integer txtL(2)
      real*8 width
      character*20 pname
 
      if(iargc().eq.0) then
         write(*,*)'Correct usage: ./main  <file with parameters>'
         write(*,*)'Example: ./main  data1.par' 
         stop
      endif

      call getarg(1,argv)                               
       
      err=readVar(argv)

      if(err.eq.-1) then
          write(*,*) 'Can not open the file'
          stop
      endif

      if(err.gt.0) then
          write(*,*)'Wrong file contents at line ',err
          stop
      endif

C  To print list od variables uncomment   
C      CALL printVar(6)

      err=sortOddParticles(cdmName)
      if(err.ne.0) then
        write(*,*) 'Can t calculate ', cdmName
        STOP
      endif 

      PDGcdm=qNumbers(cdmName,spin2, charge3,  cdim)
      write(*,fmt="('Dark matter candidate is ', A5,'with spin=',
     >I1,'/2 mass=',1PE10.2)") cdmName, spin2,Mcdm

      if(charge3.ne.0) then 
         write(*,fmt="('DM has electric charge',I1,'/3')") charge3
         STOP
      endif

      if(cdim.ne.1) then 
         write(*,*) 'DM has color charge'
         STOP
      endif


#ifdef MASSES_INFO  
      CALL printHiggs(6)
      write(*,*) '=== MASSES OF HIGGS AND ODD SECTOR: ===' 
      CALL printMasses(6,1)
#endif

#ifdef CONSTRAINTS

#endif

#ifdef HIGGSBOUNDS
      if(access(HIGGSBOUNDS//"/HiggsBounds",'x') .ne.0) call system('cd '//HIGGSBOUNDS//'; ./configure; make ')
      err=HBblocks("HB.in")
      call system(HIGGSBOUNDS// '/HiggsBounds  LandH SLHA 1 0 HB.in HB.out > hb.stdout')
      err=slhaRead("HB.out",1+4)
      write(*,fmt='("HB result=",1PE10.1,3x,"obsratio=", E10.1)')
     .slhaVal2("HiggsBoundsResults",0.d0,1,2), slhaVal2("HiggsBoundsResults",0.,1,3)
      if(0.eq.slhaSTRFormat("HiggsBoundsResults","1 5 ||%[^|]||",buff)) write(*,*) "Channel: ", buff
#endif

#ifdef HIGGSSIGNALS
#define DataSet " latestresults "
C#define Method  " peak " 
C#define  Method " mass "
#define  Method " both "
C* PDF 1 - box ; 2 = Gaussian; 3 box+Gaussia
#define PDF  " 2 " 
#define dMh " 2 " 
      write(*,*) 'HiggsSignals:'       
      if(access(HIGGSSIGNALS// '/HiggsSignals','x').ne.0 ) call system( 'cd ' //HIGGSSIGNALS//'; ./configure; make')
      call system('rm -f HS.in HS.out')
      err=HBblocks("HS.in")
      call system(HIGGSSIGNALS//'/HiggsSignals'//DataSet//Method//PDF//'SLHA 1 0 HS.in > hs.stdout');
      call system('grep -A 10000  HiggsSignalsResults HS.in > HS.out')
      err=slhaRead('HS.out',1+4)
      write(*,fmt='("  Number of observables ",F5.0)') slhaVal1("HiggsSignalsResults",0.d0,7)
      write(*,fmt='("  total chi^2=", 1PE10.1)') slhaVal1("HiggsSignalsResults",0.d0,12)
      write(*,fmt='("  p-value =", 1PE10.1)') slhaVal1('HiggsSignalsResults',0.d0,13)
#undef dMh
#undef PDF
#undef Method
#undef DataSet
#endif


#ifdef OMEGA
C to exclude processes with virtual W/Z in DM   annihilation       
      call  setVVdecay(0,0)                                 

C to include processes with virtual W/Z  also  in co-annihilation  
C     call setVVdecay(2,2)           

      write(*,*)
      write(*,*) '==== Calculation of relic density ====='
      write(*,*) 'Dark Matter candidate is ', cdmName                                                                                                             
      Omega=darkOmega(Xf,fast,Beps)
      write(*,'(a,3x,1PE10.2)') 'Omega=',Omega

      OmegaAppr=printChannels(Xf,cut,Beps,1,6)
C to restore default VZdecay and VWdecay 
       call  setVVdecay(1,1)
#endif


#ifdef  INDIRECT_DETECTION
      write(*,*) 
      write(*,*) '=============== Indirect detection ==========='

C** Emin - energy cut in GeV
      Emin=1     
      Etest=Mcdm/2
C** fi angle of sight [rad], dfi - half of cone angle    
      fi=0.1
      dfi=0.05
      sigmaV=calcSpectrum(4,SpA,SpE,SpP,NULL,NULL,NULL,err)    
C  Returns  sigma*v in cm^3/sec.
C  SpX - calculated spectra of annihilation.
C       Use SpectdNdE(E, SpX) to calculate energy distribution in  1/GeV units.
C       
C       First parameter 1-includes W/Z polarization
C                       2-includes gammas for 2->2+gamma
C                       4-print cross sections             


      write(*,*)
      write(*,fmt="('sigma*v=',1PE8.2,'[cm^3/s]')") sigmaV

      write(*,*) 'Photons'
      call gammaFluxTab(fi,dfi, sigmaV, SpA, FluxA)
      write(*,fmt="('Phonton flux  for angle of sight ',
     < F4.2,'[rad]')") fi
      write(*,fmt="('and spherical region described by cone with angle '
     >,F4.2,'[rad]' )") 2*dfi
#ifdef SHOWPLOTS      
      write(mess,fmt="('gamma flux[cm^2 s GeV]^{-1} at fi=',F4.2,
     > ' from cone with angle=',F4.2,'[rad]')")fi,dfi*2  
      err= displaySpectrum(FluxA,mess,Emin,Mcdm,1)
#endif
      write(*,fmt="('Photon flux =',1PE8.2,
     >'[cm^2 s GeV]^{-1} for E=',0PF5.1,'[GeV]')")
     > SpectdNdE(Etest, FluxA), Etest
 
      write(*,*) 'Positrons'
      call  posiFluxTab(Emin, sigmaV, SpE,  FluxE)

#ifdef SHOWPLOTS      
      err=displaySpectrum(FluxE,'positron flux [cm^2 s sr GeV]^{-1}',
     >Emin,Mcdm,1)
#endif
      write(*,fmt="('Positron flux  = ',1PE8.2,
     >'[cm^2 sr s GeV]^{-1} for E=',0PF6.1,'[GeV]')")
     > SpectdNdE(Etest, FluxE),  Etest     

      write(*,*) 'antiProtons'
      call  pbarFluxTab(Emin, sigmaV, SpP,  FluxP)

#ifdef SHOWPLOTS      
      err=displaySpectrum(FluxP,'a-proton flux [cm^2 s sr GeV]^{-1}',
     >Emin,Mcdm,1)
#endif

      write(*,fmt="('Anti proton flux  = ',1PE8.2,
     >'[cm^2 sr s GeV]^{-1} for E=',0PF6.1,'[GeV]')")
     > SpectdNdE(Etest, FluxP),  Etest     

#endif


#ifdef RESET_FORMFACTORS

C   The user has approach to form factors  which specifies quark contents 
C   of  proton and nucleon via global parametes like
C      <Type>FF<Nucleon><q>
C   where <Type> can be "Scalar", "pVector", and "Sigma"; 
C         <Nucleon>     "P" or "N" for proton and neutron
C         <q>            "d", "u","s"
C
C   calcScalarFF( Mu/Md, Ms/Md, sigmaPiN[MeV], sigma0[MeV])  
C   calculates and rewrites Scalar form factors
      write(*,*) 
      write(*,*) '===== RESET_FORMFACTOR ============'
      write(*,fmt='(" protonFF (default)  d, u, s ", 1PE10.2,E10.2,E10.2)') ScalarFFPd,ScalarFFPu,ScalarFFPs
      write(*,fmt='(" neutronFF(default)  d, u, s ", 1PE10.2,E10.2,E10.2)') ScalarFFNd,ScalarFFNu,ScalarFFNs

      write(*,*) '===========Redefinition of form factors========='  
      call calcScalarFF(0.553d0,18.9d0,70.d0,35.d0)

C    To restore default form factors of  version 2  call 
C      call  calcScalarQuarkFF(0.553d0,18.9d0,55.d0,243.5d0)


      write(*,fmt='(" protonFF (new)  d, u, s ", 1PE10.2,E10.2,E10.2)') ScalarFFPd,ScalarFFPu,ScalarFFPs
      write(*,fmt='(" neutronFF(new)  d, u, s ", 1PE10.2,E10.2,E10.2)') ScalarFFNd,ScalarFFNu,ScalarFFNs

#endif


#ifdef  CDM_NUCLEON
      write(*,*)
      write(*,*) '==== Calculation of CDM-nucleons amplitudes  ====='   

      err=nucleonAmplitudes(NoLoop, pA0,pA5,nA0,nA5)
      write(*,*) cdmName,  ' -nucleon amplitudes:'
      write(*,fmt="(' proton:  SI=',1PE9.2,' SD=',E9.2)")pA0(1),pA5(1)
      write(*,fmt="(' neutron: SI=',1PE9.2,' SD=',E9.2)")nA0(1),nA5(1) 

C  nucleon mass
      Nmass=0.939d0

      SCcoeff=4/M_PI*3.8937966E8
     >*(Nmass*Mcdm/(Nmass+ Mcdm))**2
      write(*,*) cdmName , '-nucleon cross sections[pb]:'
    
      write(*,fmt="(' proton  SI=',1PE9.2,' SD=',E9.2)")
     >                           SCcoeff*pA0(1)**2,3*SCcoeff*pA5(1)**2
      write(*,fmt="(' neutron SI=',1PE9.2,' SD=',E9.2)")
     >                           SCcoeff*nA0(1)**2,3*SCcoeff*nA5(1)**2

      write(*,*) 'anti-',cdmName,  ' -nucleon amplitudes:'
      write(*,fmt="(' proton:  SI=',1PE9.2,' SD=',E9.2)")pA0(2),pA5(2)
      write(*,fmt="(' neutron: SI=',1PE9.2,' SD=',E9.2)")nA0(2),nA5(2) 


      write(*,*) 'anti-',cdmName , '- nucleon cross sections:'
    
      write(*,fmt="(' proton  SI=',1PE9.2,' SD=',E9.2)")
     >                           SCcoeff*pA0(2)**2,3*SCcoeff*pA5(2)**2
      write(*,fmt="(' neutron SI=',1PE9.2,' SD=',E9.2)")
     >                           SCcoeff*nA0(2)**2,3*SCcoeff*nA5(2)**2

#endif


#ifdef CDM_NUCLEUS
      write(*,*)

      write(*,*) '========= Direct Detection ==============='

      write(*,*) 'CDM-73Ge:'
      nEvents=nucleusRecoil(Maxwell,73,Z_Ge,J_Ge73,
     >SxxGe73,NoLoop,dNdE)
      write(*,fmt="('73Ge: Total number of events=',1PE10.2,
     >' /day/kg')") nEvents
      write(*,fmt="(' Number of events in 10 - 50 KeV region=',
     >1PE10.2,' /day/kg')") cutRecoilResult(dNdE,10.D0,50.D0)

#ifdef SHOWPLOTS
       err=displayRecoilPlot(dNdE,
     >               'Distribution of recoil energy of 73Ge',0,199)
#endif
       write(*,*) 
       write(*,*) 'CDM-131Xe:' 
       nEvents=nucleusRecoil(Maxwell,131,Z_Xe,J_Xe131,
     > SxxXe131,NoLoop,dNdE)
      write(*,fmt="('131Xe: Total number of events=',1PE10.2,
     >' /day/kg')") nEvents
      write(*,fmt="(' Number of events in 10 - 50 KeV region=',
     >1PE10.2,' /day/kg')") cutRecoilResult(dNdE,10.D0,50.D0)

#ifdef SHOWPLOTS
      err=displayRecoilPlot(dNdE,
     >    'Distribution of recoil energy of 131Xe',0,199);
#endif

#endif

#ifdef NEUTRINO
      forSun=1 
      Emin=1  
C rhoDet - detector density
      rhoDet=1     
      write(*,*)  
      write(*,fmt='("===============Neutrino Telescope=======  for  ",$)')

      if(forSun.ne.0) then
         write(*,*) "Sun" 
      else 
         write(*,*) "Earth"
      endif 

      err=neutrinoFlux(Maxwell,forSun, FluxN,FluxNb)
#ifdef SHOWPLOTS
      err=displaySpectrum(FluxN, "nu flux from Sun [1/Year/km^2/GeV]",Emin,Mcdm,1)
      err=displaySpectrum(FluxNb,"nu-bar from Sun  [1/Year/km^2/GeV]",Emin,Mcdm,1)
#endif
      call spectrInfo(Emin/Mcdm,FluxN, IntFlux,NULL)
      write(*,fmt='("E>", 1PE8.2," GeV neutrino flux =", E8.2," [1/Year/km^2]")') Emin,IntFlux
      call  spectrInfo(Emin/Mcdm,FluxNb, IntFlux,NULL)
      write(*,fmt='("E>", 1PE8.2," GeV anti-neutrino flux =", E8.2," [1/Year/km^2]")') Emin,IntFlux
  
C ======  Upward events
  
      call muonUpward(FluxN,FluxNb,FluxMu)
#ifdef SHOWPLOTS  
      err= displaySpectrum(FluxMu,"Upward muons[1/Year/km^2/GeV]",Emin,Mcdm/2,1)
#endif

      call spectrInfo(Emin/Mcdm,FluxMu, IntFlux,NULL)
      write(*,fmt='("E>", 1PE8.2," GeV Upward muon flux =", E8.2," [1/Year/km^2]")') Emin,IntFlux
  
C ======  Contained events
      call muonContained(FluxN,FluxNb,rhoDet, FluxMu)
#ifdef SHOWPLOTS  
      err=displaySpectrum(FluxMu,"Contained  muons[1/Year/km^3/GeV]",Emin,Mcdm,1)
#endif
      call spectrInfo(Emin/Mcdm,FluxMu, IntFlux,NULL)
      write(*,fmt='("E>", 1PE8.2," Contained muon flux =", E8.2," [1/Year/km^3]")') Emin,IntFlux
#endif

#ifdef DECAYS

      write(*,*)
      write(*,*)'======Calculations of widths======='

      write(*,*)
      pname = 'h'
      width=pWidth(pname,txtL)
      write(*,fmt="(A8, ':  total width ',1PE8.2,'[GeV]')") pname, width
      write(*,*)'and Branchings:'
      call printTxtList(txtL,6)  

      write(*,*)
      pname = '~H+'
      width=pWidth(pname,txtL)
      write(*,fmt="(A8,':  total width ',1PE8.2,'[GeV]')")pname, width
      write(*,*)'and Branchings:'
      call printTxtList(txtL,6)  

      write(*,*)
      pname = '~H3'
      width=pWidth(pname,txtL)
      write(*,fmt="(A8,':  total width ',1PE8.2,'[GeV]')")pname, width
      write(*,*)'and Branchings:'
      call printTxtList(txtL,6)  

#endif
#ifdef CLEAN
      call system("rm -f HB.in HB.out HS.in HS.out hb.stdout hs.stdout debug_channels.txt debug_predratio.txt Key.dat")
#endif

      call killPlots
      END
