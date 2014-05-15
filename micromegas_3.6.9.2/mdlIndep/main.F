
C====== Modules =================
C   Keys to switch on
C   various modules of micrOMEGAs
C================================
   
#define INDIRECT_DETECTION
C   Compute spectra of gamma/positron/antiprotons/neutrinos for DM annihilation; 
C   Input parameters are  Mcdm,  annihilation cross section vcs and channel fractions.
C   Integrate gamma signal over DM galactic squared density for given line of sight; 
C   Calculate galactic propagation of positrons and antiprotons.      
#define CDM_NUCLEUS
C   Calculate number of events for 1kg*day and recoil energy distibution for various nuclei.
C   Input parameters are  Mcdm, spin dependent and spin independent cross sections on nucleons
C================ end of Modules  ======

#define NEUTRINO


C#define SHOWPLOTS
C     Display  graphical plots on the screen

#define CLEAN to crean intermediate files

C===== End of DEFINE  settings ===== 


      IMPLICIT NONE
#include "../sources/micromegas.fh"
#include "lib/pmodel.fh"

C  =================      INPUT PARAMETERS  =============
C number of annihilation channels
#define  nCH 5   
C PDG of annihilation channels
      integer    pdgCH(nCH) /   6,   5,  15,  23,   24/   
C  fraction of annihilation channels    
      real*8     fracCH(nCH) /0.59,0.12,0.20,0.03, 0.06/
C DM annihilation cross section times  velocity
      real *8 vcs/0.310/ 
      real*8 csSIp/5.489E-09/,csSIn/5.758E-09/,csSDp/5.807E-05/,csSDn/-4.503E-05/

C ============ Local Variables 
C  for Indirect Detection    
      real*8  Emin,sigmaV,fi,dfi,Etest
      real*8 SpA(NZ),Spe(NZ),SpP(NZ) 
      real*8 FluxA(NZ),FluxE(NZ),FluxP(NZ)
      real*8 buff(NZ)
      real*8 SMmev/320/ 
      integer i,n,err
C  for CDM_NUCLEUS
      real*8  dNdE(200),rho,nEvents

C  for NEUTRINO 
      integer forSun
      real*8 IntFlux, rhoDet
      real*8 nu(NZ),nu_bar(NZ),FluxN(NZ),FluxNb(NZ),FluxMu(NZ) 
      real*8 Prop,cRate,R       

      Mcdm=189.0

#ifdef  INDIRECT_DETECTION
      write(*,*) 
      write(*,*) '=============== Indirect detection ==========='
C** Emin - energy cut in GeV
      Emin=1     
      Etest=Mcdm/2
C** fi angle of sight [rad], dfi - half of cone angle    
      fi=0.
      dfi=3.141592654/180.
      write(*,*)
      sigmaV=vcs*2.9979E-26;
      write(*,fmt="('sigma*v=',1PE8.2,'[cm^3/s]')") sigmaV 
      DO i=1,NZ
        SpA(i)=0
        SpE(i)=0
        SpP(i)=0
      ENDDO

      DO  n=1,nCH
        if(fracCH(n).gt.0) then
          err=basicSpectra(pdgCH(n),0, buff)
          DO i=1,NZ
            SpA(i)=SpA(i)+buff(i)*fracCH(n)
          ENDDO
          err= basicSpectra(pdgCH(n),1, buff)
          DO  i=1,NZ
            SpE(i)=SpE(i)+buff(i)*fracCH(n)
          ENDDO  
          err= basicSpectra(pdgCH(n),2, buff)
          DO i= 1,NZ
            SpP(i)= SpP(i)+buff(i)*fracCH(n) 
          ENDDO
        ENDIF
      ENDDO
      write(*,*) 'Photons'
      call gammaFluxTab(fi,dfi, sigmaV, SpA, FluxA)
      write(*,fmt="('Phonton flux  for angle of sight ',F4.2,'[rad]')") fi
      write(*,fmt="('and spherical region described by cone with angle ',F4.2,'[rad]' )") 2*dfi
#ifdef SHOWPLOTS      
      write(mess,fmt="('gamma flux[cm^2 s GeV]^{-1} at fi=',F4.2,' from cone with angle=',F4.2,'[rad]')")fi,dfi*2  
      err= displaySpectrum(FluxA,mess,Emin,Mcdm,1)
#endif
      write(*,fmt="('Photon flux =',1PE8.2,'[cm^2 s GeV]^{-1} for E=',0PF5.1,'[GeV]')")
     > SpectdNdE(Etest, FluxA), Etest
 
      write(*,*) 'Positrons'
      call  posiFluxTab(Emin, sigmaV, SpE,  FluxE)
      if(SMmev.gt.0)  call solarModulation(SMmev,5D-4,FluxE,FluxE)
#ifdef SHOWPLOTS      
      err=displaySpectrum(FluxE,'positron flux [cm^2 s sr GeV]^{-1}',
     >Emin,Mcdm,1)
#endif
      write(*,fmt="('Positron flux  = ',1PE8.2,
     >'[cm^2 sr s GeV]^{-1} for E=',0PF6.1,'[GeV]')")
     > SpectdNdE(Etest, FluxE),  Etest     

      write(*,*) 'antiProtons'
      call  pbarFluxTab(Emin, sigmaV, SpP,  FluxP)
      if(SMmev.gt.0)  call solarModulation(SMmev,1.d0,FluxP,FluxP)
#ifdef SHOWPLOTS      
      err=displaySpectrum(FluxP,'a-proton flux [cm^2 s sr GeV]^{-1}',
     >Emin,Mcdm,1)
#endif

      write(*,fmt="('Anti proton flux  = ',1PE8.2,
     >'[cm^2 sr s GeV]^{-1} for E=',0PF6.1,'[GeV]')")
     > SpectdNdE(Etest, FluxP),  Etest     

#endif



#ifdef CDM_NUCLEUS
      write(*,*)

      write(*,*) '========= Direct Detection ==============='

      write(*,*) 'CDM-73Ge:'
      nEvents=nucleusRecoilAux(Maxwell,73,Z_Ge,J_Ge73,
     >SxxGe73,csSIp,csSIn, csSDp,csSDn,dNdE)
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
       nEvents=nucleusRecoilAux(Maxwell,131,Z_Xe,J_Xe131,
     > SxxXe131,csSIp,csSIn, csSDp,csSDn,dNdE)
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

C      err=neutrinoFlux(Maxwell,forSun, FluxN,FluxNb)
       Crate=captureAux(forSun,Maxwell, csSIp,csSIn,csSDp,csSDn)
 
       DO i=1,NZ
         nu(i)=0
         nu_bar(i)=0
       ENDDO    

       DO n=1,nCH
         if(fracCH(n).gt.0) then 
           err=  basicNuSpectra(forSun, pdgCH(n), 1, buff) 
           DO i=1,NZ
              nu(i) = nu(i)+buff(i)*fracCH(n)
           ENDDO
           err=basicNuSpectra(forSun, pdgCH(n),-1, buff)
           DO  i=1,NZ 
             nu_bar(i)=nu_bar(i)+buff(i)*fracCH(n)
           ENDDO
         endif
       ENDDO  


C distance in km
       if(forSun.ne.0) then 
          R=150E6 
       else  
          R=6378.1
       endif
C flux normalization in  in Year*km^2
       Prop=31556925.2/(4*M_PI*R*R)    
       DO i=1, NZ 
         FluxN(i)=nu(i)*0.5*Crate*Prop   
         FluxNb(i)=nu_bar(i)*0.5*Crate*Prop 
       ENDDO

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

      call pdg2name(23,pname)
      if(pname .ne. ' ') then
         width=pWidth(pname,txtL)
         write(*,fmt="(A5,':  total width ',1PE8.2,'[GeV]')") pname,width
         write(*,*)'and Branchings:'
         call printTxtList(txtL,6)
      endif 

      call pdg2name(25,pname)
      if(pname .ne. ' ') then
         width=pWidth(pname,txtL)
         write(*,fmt="(A5,' :  total width ',1PE8.2,'[GeV]')") pname,width
         write(*,*)'and Branchings:'
         call printTxtList(txtL,6)
      endif 

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

#endif

      call killPlots
      END
