C=========   NMSSM scenario  ==========
C  One can define SUGRA or EWSB 
C  Default main.F reads SLHA file.
C=====================================

C#define SUGRA
#define EWSB

C====== Modules =================
C   Keys to switch on
C   various modules of micrOMEGAs
C================================
   
#define MASSES_INFO
C    Display information about SUSY and Higgs masses
#define CONSTRAINTS
C     Display  deltarho, B_>sgamma, Bs->mumu, gmuon and
C     check LEP mass limits
 
C#define HIGGSBOUNDS  '../Packages/HiggsBounds-4.1.0'
C#define HIGGSSIGNALS '../Packages/HiggsSignals-1.2.0'

#define OMEGA
C    Calculate relic density and display contribution of
C                            individual channels 
#define INDIRECT_DETECTION
C     Compute spectra of gamma/positron/neutrinos
C     for DM annihilation; calculate <sigma*v> and
C     integrate gamma signal over DM galactic squared
C     density for given line of sight.

C#define LoopGAMMA
C     Calculate discrete  photon spectrum caused by annihilation of
C      neutralinos into two photons and Z-photon


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

#define NEUTRINO for neutrini telescope

#define DECAYS
C     Calculates particle widths and branchings

C#define CROSS_SECTIONS
C     Calculate cross sections and widths for
C      reactions specified by the user
 
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

C! local  Parameters   
      INTEGER fast/1/
      REAL*8  Beps/1.E-5/,cut/0.01/
      REAL*8  m0,mhf,a0,tb,sgn
      REAL*8 m32

C Variables      
      INTEGER err
      REAL*8 Omega,OmegaAppr,Xf,val
      REAL*8  gMG1,gMG2,gMG3,gAl,gAt,gAb,gMHu,gMHd,
     > gMl2,gMl3,gMr2,gMr3,gMq2,gMq3,gMu2,gMu3,gMd2,gMd3,
     >Lambda,aKappa,aLambda
      INTEGER  narg,nw 
      character *100 argv, buff
      character*10 cdmName
      integer spin2,charge3,cdim,PDGcdm
      real*8 Pcm,cosmin, cosmax,rw
      integer cc(2)
      integer ntot,l, nin,nout
      character*6 name(4)
      real*8 mass(4)
      real*8 Etest,cs
      integer outP
      real*8 Emin,NGtot,Xtot,sigmaV,v/0.001/,tab(NZ),fi,dfi
      character*100 mess

C  for CONSTRAINTS    
      real*8 constr0,constrM, constrP

C  for INDIRECT 
      real*8 N_tot, SpA(NZ),SpE(NZ), SpP(NZ) 
      real*8 FluxA(NZ),FluxE(NZ), FluxP(NZ)
C for LoopGAMMA
      real*8 vcs_gg,vcs_zg
C for CDM_NUCLEON
      real*8  pA0(2),pA5(2),nA0(2),nA5(2),Nmass,SCcoeff,
     >ffS0P(3),ffS0N(3),ffbuff(3)

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

C   to save the SLHA  input/output file  uncomment  


#ifdef SUGRA
      write(*,*)
      write(*,*) '========= mSUGRA NMSSM scenario ========'
      narg = iargc()
      if(narg.lt.6) then
        write(*,*) ' One has to launch this program with parameters.'
        write(*,*) ' Six  of them are obligatory:'
        write(*,*) '    m0      common scalar mass at GUT scale'
        write(*,*) '    mhf     common gaugino mass at GUT scale'
        write(*,*) '    a0      trilinear soft breaking at GUT scale'
        write(*,*) '    tb      tan(beta)'
        write(*,*) '   Lambda   Lambda parameter'
        write(*,*) '   aKappa   aKappa parameter'
        write(*,*) ' Auxiliary parameters are'
        write(*,*) '    sgn  1(default)/-1, sign of Higgsino mass term'
        write(*,*) '    Mtp     top quark pole mass'
        write(*,*) '    MbMb    Mb(Mb) scale independent b-quark mass'
        write(*,*) '    alfSMZ  strong coupling at MZ'
        write(*,*) ' ================================================'
        write(*,*) 'Example: ./main 320 600 -1300 2 0.5 -1400'
        STOP
      else
        call getarg(1, argv)
        read(argv,*,ERR=100) m0
        call getarg(2, argv)
        read(argv,*,ERR=100) mhf
        call getarg(3, argv)
        read(argv,*,ERR=100) a0
        call getarg(4, argv)
        read(argv,*,ERR=100) tb
        call getarg(5, argv)
        read(argv,*,ERR=100) Lambda
        call getarg(6, argv)
        read(argv,*,ERR=100) aKappa

        if(narg.ge.7) then
          call getarg(7, argv)
          read(argv,*,ERR=100) sgn
        else 
          sgn=1
        endif 
        if(narg.ge.8) then
          call getarg(8, argv)
          read(argv,*,ERR=100) aLambda
        else 
          aLambda=a0 
        endif 

        if(narg.ge.9) then
          call getarg(9, argv)
          read(argv,*,ERR=100) val
          call assignValW('Mtp',val)
        endif 
        if(narg.ge.10) then
          call getarg(10, argv)
          read(argv,*,ERR=100) val
          call assignValW('MbMb',val)
        endif 
        if(narg.ge.11) then
          call getarg(11, argv)
          read(argv,*,ERR=100) val
          call assignValW('alfSMZ',val)
        endif 
      endif
  
      err=nmssmSUGRA(m0,mhf,a0,tb,sgn,Lambda,aLambda,aKappa)
      if(err.ne.0) goto 100

#elif defined(EWSB)

      write(*,*)
      write(*,*) '========= EWSB scale input NMSSM========='
      if(iargc().eq.0) then
         write(*,*)'Correct usage: ./main  <file with NMSSM parameters>'
         write(*,*)'Example: ./main  data1.par' 
         stop
      endif

      call getarg(1,argv)                               

      err=readVarNMSSM(argv)

      if(err.eq.-1) then
          write(*,*) 'Can not open the file'
          stop
      endif

      if(err.gt.0) then
          write(*,*)'Wrong file contents at line ',err
          stop
      endif

      err=nmssmEWSB(0)
      if(err.lt.0) goto 100
#else
      write(*,*)
      write(*,*) '========= SLHA file input ========='
      narg=iargc()
      if(narg.eq.0) then 
        write(*,*) 
     >  'The program needs one argument:the name of SLHA input file.'
        write(*,*)  'Example: ./main spectr.dat'
        stop
      endif   
      call getarg(1,argv)
      write(*,*) 'Initial file  ', argv
      err=lesHinput(argv)
      write(*,*) 'err=', err
      if(err.ne.0) stop

#endif

      goto 101
100   STOP '!Wrong input' 
101   continue
C  To print list od variables uncomment   
C      CALL printVar(6)

C    to save results in SLHA format  use 
C     Call writeLesH('slha.out') 


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
      CALL printHiggs(6)
      CALL printMasses(6,1)
#endif

#ifdef CONSTRAINTS
      write(*,*)

      IF(err.ne.0) THEN
         WRITE(*,*)' No RGE solution'
         WRITE(*,*)' Error message form NMHTools'
         nw=NMHwarn(6)
         STOP
      ELSE 
        write(*,*) 'Warnings from NMSSMTools:'
        nw=NMHwarn(6)
        if(nw.eq.0) write(*,*) '  .......  none'
      ENDIF
      write(*,*)
      write(*,*) '==== Physical Constraints: ====='

      constr0=bsgnlo(constrM,constrP)
      write(*,fmt="('B->s,gamma = ',1PE9.2,'(',E9.2,',',E9.2,')')") constr0,constrM, constrP 

      constr0= bsmumu(constrM,constrP)
      write(*,fmt="('Bs->mu,mu  = ',1PE9.2,'(',E9.2,',',E9.2,')')") constr0,constrM, constrP  
  
      constr0=btaunu(constrM,constrP)
      write(*,fmt="('B+->tau+,nu= ',1PE9.2,'(',E9.2,',',E9.2,')')") constr0,constrM, constrP 

  
      constr0=deltamd(constrM,constrP)
      write(*,fmt="('deltaMd    = ',1PE9.2,'(',E9.2,',',E9.2,')')")  constr0,constrM, constrP 

      constr0=deltams(constrM,constrP)
      write(*,fmt="('deltaMs    = ',1PE9.2,'(',E9.2,',',E9.2,')')") constr0,constrM, constrP 

      constr0=gmuon(constrM,constrP)
      write(*,fmt="('gmuon/BSM  = ',1PE9.2,'(',E9.2,',',E9.2,')')") constr0,constrM, constrP 
#endif

#ifdef HIGGSBOUNDS
      if(access(HIGGSBOUNDS//"/HiggsBounds",'x') .ne.0) call system('cd '//HIGGSBOUNDS//'; ./configure; make ')
      call system('cat spectr decay > HB.in')
      err=HBblocks('HB.in')
      call system( HIGGSBOUNDS//'/HiggsBounds  LandH SLHA 5 1 HB.in HB.out >hb.stdout')
      err=slhaRead("HB.out",1+4)
      write(*,fmt='("HB result=",1PE10.1,"obsratio=", E10.1)')
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
      err=slhaWrite("HS.in")
      err=HBblocks("HS.in")
      call system('echo "BLOCK DMASS" >> HS.in')
      call system('echo "  25 '//dMh// '" >> HS.in')
      call system(HIGGSSIGNALS//'/HiggsSignals'//DataSet//Method//PDF//'SLHA 5 1 HS.in > hs.stdout');
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
      if(cdmName.ne.'~o1') then
          write(*,*) '~o1 is not LSP'
      else
          call o1Contents(6)
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
      call gammaFluxTab(fi,dfi, sigmaV, SpA,  FluxA)
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
      call  posiFluxTab(Emin, sigmaV, SpE, FluxE)

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

#ifdef LoopGAMMA
        
      if(loopGamma(vcs_zg, vcs_gg).eq.0) then 
         write(*,*) '=========== Gamma ray lines =============:'
         write(*,*)    
         write(*,fmt="('E=',F6.1,'[GeV] vcs(Z,A)=', 1PE8.2,
     >        '[cm^2/s] flux=',E8.2,'[cm^2 s]^{-1}')") 
     >   Mcdm-91.19*91.19/4/Mcdm, vcs_zg, 
     >   gammaFlux(fi,dfi,vcs_zg)

         write(*,fmt="('E=',F6.1,'[GeV] vcs(A,A)=', 1PE8.2,
     >        '[cm^2/s] flux=',E8.2,'[cm^2 s]^{-1}')") 
     >   Mcdm, vcs_gg, 
     >   2*gammaFlux(fi,dfi,vcs_gg)

      endif
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
      write(*,*)'CDM-nucleon amplitudes:'
      write(*,fmt="(' proton:  SI=',1PE9.2,' SD=',E9.2)")pA0(1),pA5(1)
      write(*,fmt="(' neutron: SI=',1PE9.2,' SD=',E9.2)")nA0(1),nA5(1) 

C  nucleon mass
      Nmass=0.939d0

      SCcoeff=4/M_PI*3.8937966E8
     >*(Nmass*Mcdm/(Nmass+ Mcdm))**2
      write(*,*) 'CDM-nucleon cross sections[pb]:'

      write(*,fmt="(' proton  SI=',1PE9.2,' SD=',E9.2)")
     >                           SCcoeff*pA0(1)**2,3*SCcoeff*pA5(1)**2
      write(*,fmt="(' neutron SI=',1PE9.2,' SD=',E9.2)")
     >                           SCcoeff*nA0(1)**2,3*SCcoeff*nA5(1)**2
#endif


#ifdef CDM_NUCLEUS
      write(*,*)
      write(*,*)'=========  Direct Detection ==========='

      write(*,*) 'Ge:'

      nEvents=nucleusRecoil(Maxwell,73,Z_Ge,J_Ge73,
     >S00Ge73,S01Ge73,S11Ge73,FeScLoop,dNdE)
C  See '../sources/micromegas.h' for description of arguments
C  Instead of Maxwell  one can use 'fDvDelta' Delta-function
C  velocity distribution.
      write(*,fmt="('Total number of events=',1PE9.2,' /day/kg')") 
     >nEvents

      write(*,fmt=
     >"(' Number of events in 10 - 50 KeV region=',1PE9.2,' /day/kg')") 
     >cutRecoilResult(dNdE,10.d0,50.d0) 


#ifdef SHOWPLOTS
       err=displayRecoilPlot(dNdE,
     >           'Distribution of recoil energy of 73Ge',0.d0,199.d0)
#endif
       write(*,*)
       write(*,*) '131Xe:'  
       nEvents=nucleusRecoil(Maxwell,131,Z_Xe,J_Xe131,
     > S00Xe131,S01Xe131,S11Xe131,FeScLoop,dNdE)

      write(*,fmt=
     >"('Total number of events=',1PE9.2,' /day/kg')")nEvents
      write(*,fmt="(' Number of events in 10 - 50 KeV region=',
     >   1PE9.2, ' /day/kg')")cutRecoilResult(dNdE,1.d1,5.d1)


#ifdef SHOWPLOTS
      err=displayRecoilPlot(dNdE,
     >    'Distribution of recoil energy of 131Xe',0.d0,199.d0)
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
      pname = 'h1'
      width=pWidth(pname,txtL)
      write(*,fmt="(A8, ':  total width ',1PE8.2,'[GeV]')") pname, width
      write(*,*)'and Branchings:'
      call printTxtList(txtL,6)  

      write(*,*)
      pname = '~o2'
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

#ifdef CLEAN 
      call system("rm -f inp spectr decay omega nngg.out  Key.dat")
      call system("rm -f debug_channels.txt debug_predratio.txt HB.in HB.out HS.in HS.out hb.stdout hs.stdout")

#endif 


        call killPlots
      END
