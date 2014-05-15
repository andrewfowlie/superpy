
C====== Modules =================
C   Keys to switch on
C   various modules of micrOMEGAs
C================================
   
#define MASSES_INFO
C    Display masses of odd sector
#define CONSTRAINTS
C    Not implemented
#define OMEGA
C    Calculate relic density and display contribution of
C                            individual channels 
#define INDIRECT_DETECTION
C     Compute spectra of gamma/positron/neutrinos
C     for DM annihilation; calculate <sigma*v> and
C     integrate gamma signal over DM galactic squared
C     density for given line of sight.
#define RESET_FORMFACTORS
C     Modify default nucleus form factors, 
C     DM velocity distribution,
C     A-dependence of Fermi-dencity.
#define CDM_NUCLEON 
C     Calculate amplitudes and cross-sections for
C     CDM-mucleon collisions 
#define CDM_NUCLEUS
C     Calculate number of events for 1kg*day
C     and recoil energy distibution for various nuclei

#define NEUTRINO for neutrino telescope

#define DECAYS
    /* Calculation of cross sections and branchings */

C#define CROSS_SECTIONS
C     Calculate cross sections and widths for
C      reactions specified by the user
 
C================ end of Modules  ======

C===== Options ========

C#define SHOWPLOTS
C     Display  graphical plots on the screen

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
      REAL*8 Omega,OmegaAppr,Xf,val
      INTEGER  narg,nw 
      character *100 argv
      character*10 cdmName
      integer spin2,charge3,cdim,PDGcdm
      real*8 Pcm,cosmin, cosmax
      integer cc(2)
      integer ntot,l, nin,nout
      character*6 name(4)
      real*8 mass(4)
      real*8 cs
      real*8  Emin,Etest,sigmaV,fi,dfi
      real*8 SpA(NZ),Spe(NZ),SpP(NZ)
      real*8 FluxA(NZ),FluxE(NZ),FluxP(NZ)
      character*100 buff

C for CDM_NUCLEON
      real*8  pA0(2),pA5(2),nA0(2),nA5(2),Nmass,SCcoeff,
     >ffS0P(3),ffS0N(3),ffbuff(3)

C  for CDM_NUCLEUS
      real*8  dNdE(200),nEvents
C  for NEUTRINO
      integer forSun
      real*8 IntFlux, rhoDet   
      real*8 FluxN(NZ),FluxNb(NZ),FluxMu(NZ)
C  For DECAYS
      integer txtL(2)
      real*8 width
      character*20 pname


      write(*,*)
      write(*,*) '========= EWSB scale input For CPVMSSM ========='
      if(iargc().eq.0) then
         write(*,*)'Correct usage: ./main  <file with parameters>'
         write(*,*)'Example: ./main  data1.par' 
         stop
      endif

      call getarg(1,argv)                               
       
      err=readVarRHNM(argv)

      if(err.eq.-1) then
          write(*,*) 'Can not open the file'
          stop
      endif

      if(err.gt.0) then
          write(*,*)'Wrong file contents at line ',err
          stop
      endif

      goto 101
100   STOP 'Wrong input' 
101   continue
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
      write(*,*) '=== MASSES OF HIGGS AND SUSY PARTICLES: ===' 
      CALL printMasses(6,1)
#endif

#ifdef CONSTRAINTS
#endif


#ifdef OMEGA

C to exclude processes with virtual W/Z in DM   annihilation       
       call  setVVdecay(0,0)                                 

C to include processes with virtual W/Z  also  in co-annihilation  
C      call setVVdecay(2,2)           

      write(*,*)
      write(*,*) '==== Calculation of relic density ====='
      write(*,*) 'Dark Matter candidate is ', cdmName                                   
      if(cdmName.ne.'~n4') then
          write(*,*) '~n4 is not LSP'
      endif                                                                          
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
      sigmaV=calcSpectrum(4+1,SpA,SpE,SpP,NULL,NULL,NULL,err)    
C  Returns  sigma*v in cm^3/sec.
C  SpX - calculated spectra of annihilation.
C  Use SpectdNdE(E, SpX) to calculate energy distribution in 1/GeV units.
      write(*,*)
      write(*,fmt="('sigma*v=',1PE8.2,'[cm^3/s]')") sigmaV

      write(*,*) 'Photons'
      call gammaFluxTab(fi,dfi, sigmaV, SpA, FluxA)
      write(*,fmt="('Phonton flux  for angle of sight ',
     < F4.2,'[rad]')") fi
      write(*,fmt="('and spherical region described by cone with angle '
     >,F4.2,'[rad]' )") 2*dfi
#ifdef SHOWPLOTS      
      write(buff,fmt="('gamma flux[cm^2 s GeV]^{-1} at fi=',F4.2,
     > ' from cone with angle=',F4.2,'[rad]')")fi,dfi*2  
      err= displaySpectrum(FluxA,buff,Emin,Mcdm,1)
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
C   calcScalarQuarkFF( Mu/Md, Ms/Md, sigmaPiN[MeV], sigma0[MeV])  
C   calculates and rewrites Scalar form factors

      write(*,fmt='(" protonFF (default)  d, u, s ", 1PE10.2,E10.2,E10.2)') ScalarFFPd,ScalarFFPu,ScalarFFPs
      write(*,fmt='(" neutronFF(default)  d, u, s ", 1PE10.2,E10.2,E10.2)') ScalarFFNd,ScalarFFNu,ScalarFFNs

      write(*,*) '===========Redefinition of form factors========='
  
      call calcScalarQuarkFF(0.46d0,27.5d0,34.d0,42.d0)

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

C  nucleon mass
      Nmass=0.939d0

      SCcoeff=4/M_PI*3.8937966E8
     >*(Nmass*Mcdm/(Nmass+ Mcdm))**2
      write(*,*) 'anti-',cdmName , '- nucleon cross sections:'
    
      write(*,fmt="(' proton  SI=',1PE9.2,' SD=',E9.2)")
     >                           SCcoeff*pA0(2)**2,3*SCcoeff*pA5(2)**2
      write(*,fmt="(' neutron SI=',1PE9.2,' SD=',E9.2)")
     >                           SCcoeff*nA0(2)**2,3*SCcoeff*nA5(2)**2
#endif


#ifdef CDM_NUCLEUS
      write(*,*)
      write(*,*) '============= Direct Detection ============'

      write(*,*) '73Ge:'
      nEvents=nucleusRecoil(Maxwell,73,Z_Ge,J_Ge73,
     >SxxGe73,NoLoop,dNdE)

      write(*,fmt="('Total number of events=',1PE9.2,' /day/kg')") 
     >nEvents

      write(*,fmt=
     >"(' Number of events in 10 - 50 KeV region=',1PE9.2,' /day/kg')") 
     >cutRecoilResult(dNdE,1.D1,5.D1) 

#ifdef SHOWPLOTS
       err=displayRecoilPlot(dNdE,
     >               'Distribution of recoil energy of 73Ge',0,199)
#endif
       write(*,*) 
       write(*,*) '131Xe:' 
       nEvents=nucleusRecoil(Maxwell,131,Z_Xe,J_Xe131,
     > SxxXe131,NoLoop,dNdE)

      write(*,fmt=
     >"('Total number of events=',1PE9.2,' /day/kg')")nEvents
      write(*,fmt="(' Number of events in 10 - 50 KeV region=',
     >   1PE9.2, ' /day/kg')")cutRecoilResult(dNdE,1.D1,5.D5)
#ifdef SHOWPLOTS
      err=displayRecoilPlot(dNdE,
     >    'Distribution of recoil energy of 131Xe',0,199)
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
      pname = 'H'
      width=pWidth(pname,txtL)
      write(*,fmt="(A8, ':  total width ',1PE8.2,'[GeV]')") pname, width
      write(*,*)'and Branchings:'
      call printTxtList(txtL,6)  

      write(*,*)
      pname = 'Zp'
      width=pWidth(pname,txtL)
      write(*,fmt="(A8,':  total width ',1PE8.2,'[GeV]')")pname, width
      write(*,*)'and Branchings:'
      call printTxtList(txtL,6)  

#endif



#ifdef CROSS_SECTIONS
      write(*,*)
      write(*,*)'======Calculations of cross sections======='

      Pcm=500
      cosmin=-0.99
      cosmax=0.99
      write(*,fmt="('Process e,E->2*x at Pcm= ',F6.1,'[GeV]')") Pcm
      call setHelicities(0.5D0,-0.5D0)
      call newProcess('e%,E%->2*x',cc)
      if(cc(1).ne.0 ) then
        call procInfo1(cc,ntot,nin,nout)
        do l=1,ntot
          cs= cs22(cc,l,Pcm,cosmin,cosmax,err)
          if((err.eq.0) .and. (cs.ne.0)) then 
            call procInfo2(cc,l,name,mass)
            write(*,fmt='(A3,1x,A3,A4,A3,1x,A3,$)')
     >      name(1),name(2),' -> ', name(3),name(4)
            write(*,fmt="(1PE8.2,'[pb]')") cs
          endif 
        enddo
      endif     

#endif
        call killPlots
      END
