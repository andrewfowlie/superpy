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
#define OMEGA
C    Calculate relic density and display contribution of
C                            individual channels 
#define INDIRECT_DETECTION
C     Compute spectra of gamma/positron/neutrinos
C     for DM annihilation; calculate <sigma*v> and
C     integrate gamma signal over DM galactic squared
C     density for given line of sight.
C!#define RESET_FORMFACTORS
C     Modify default nucleus form factors, 
C     DM velocity distribution,
C     A-dependence of Fermi-dencity.
#define CDM_NUCLEON 
C     Calculate amplitudes and cross-sections for
C     CDM-mucleon collisions 
C#define CDM_NUCLEUS
C     Calculate number of events for 1kg*day
C     and recoil energy distibution for various nuclei

C#define DECAYS
C     Calculates particle widths and branchings

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
      character *100 argv
      character*10 cdmName
      integer spin2,charge3,cdim,PDGcdm
      real*8 Pcm,cosmin, cosmax,rw,width
      integer cc(2),txtL(2),dim
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
C for CDM_NUCLEON
      real*8  pA0(2),pA5(2),nA0(2),nA5(2),Nmass,SCcoeff,
     >ffS0P(3),ffS0N(3),ffbuff(3)

C  for CDM_NUCLEUS
      real*8  dNdE(200),rho,nEvents

C   to save the SLHA  input/output file  uncomment  
c       call delFiles(0)

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

      err=nmssmEWSB()
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
      err=readSLHA(argv)
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
      write(*,fmt="('B->s,gamma = ',1PE9.2,'(',E9.2,',',E9.2,')')") 
     >constr0,constrM, constrP 

      constr0= bsmumu(constrM,constrP)
      write(*,fmt="('Bs->mu,mu  = ',1PE9.2,'(',E9.2,',',E9.2,')')") 
     >constr0,constrM, constrP  
  
      constr0=btaunu(constrM,constrP)
      write(*,fmt="('B+->tau+,nu= ',1PE9.2,'(',E9.2,',',E9.2,')')") 
     >constr0,constrM, constrP 

  
      constr0=deltamd(constrM,constrP)
      write(*,fmt="('deltaMd    = ',1PE9.2,'(',E9.2,',',E9.2,')')") 
     >constr0,constrM, constrP 

      constr0=deltams(constrM,constrP)
      write(*,fmt="('deltaMs    = ',1PE9.2,'(',E9.2,',',E9.2,')')") 
     >constr0,constrM, constrP 

      constr0=gmuon(constrM,constrP)
      write(*,fmt="('gmuon/BSM  = ',1PE9.2,'(',E9.2,',',E9.2,')')") 
     >constr0,constrM, constrP 
 
#endif


#ifdef OMEGA
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

      write(*,*) 'protonFF (default)  d, u, s ', 
     >       ScalarFFPd, ScalarFFPu,ScalarFFPs
      write(*,*) 'neutronFF(default)  d, u, s ', 
     >       ScalarFFNd, ScalarFFNu,ScalarFFNs                          
      write(*,*) '===========Redefinition of form factors========='  
      call calcScalarQuarkFF(0.553d0,18.9d0,45.5d0,26.d0)
      write(*,*) 'protonFF (new)      d, u, s ', 
     >ScalarFFPd, ScalarFFPu,ScalarFFPs                           
      write(*,*) 'neutronFF(new)      d, u, s ', ScalarFFNd, 
     >ScalarFFNu,ScalarFFNs                    

C  Option to change parameters of DM velocity  distribution 

      call SetfMaxwell(220.d0,600.d0)
C    dN  ~  exp(-v^2/arg1^2)*Theta(v-arg2)  d^3v     
C    Earth velocity with respect to Galaxy defined by 'Vearth' parameter.
C    All parameters are  in [km/s] units.       
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

#ifdef DECAYS

      write(*,*)
      write(*,*)'======Calculations of widths======='

C      pname = "h";
      width=pWidth('h1',txtL,dim)
      write(*,fmt="('h1->',I1,'*x :  total width ',1PE8.2,'[GeV]')")
     >dim,width
      write(*,*)'and Branchings:'
      call printTxtList(txtL,6)

      width=pWidth('l',txtL,dim)
      write(*,fmt="('l->',I1,'*x :  total width ',1PE8.2,'[GeV]')")
     >dim,width
      write(*,*)'and Branchings:'
      call printTxtList(txtL,6)

      write(*,fmt="('Br(e,Ne,nl)=',F4.2)") findBr(txtL,"e,Ne,nl")

C   pname = "~o2";
      width=pWidth('~o2',txtL,dim)
      write(*,fmt="('~o2->',I1,'*x :  total width ',1PE8.2,'[GeV]')")
     >dim,width
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
      call newProcess('e%,E%->2*x','eE_2x',cc)
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
