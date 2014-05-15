! This file is part of HiggsBounds
!  -KW
!******************************************************
      program SM_vs_fourth_generation_model
!******************************************************
      implicit none
      integer nHzero,nHplus,HBresult,chan,ncombined,i
      real*8 obsratio,mass
      real*8 SMGammaTotal,SMGamma_hgg,Gamma_hgg,SMGamma_h,
     &                SMBR_Hgg,FGMGammaTotal

      real*8 Mh,GammaTotal,g2hjss_s,g2hjss_p,g2hjcc_s,g2hjcc_p,
     &       g2hjbb_s,g2hjbb_p,g2hjtt_s,g2hjtt_p,
     &       g2hjmumu_s,g2hjmumu_p,g2hjtautau_s,g2hjtautau_p,
     &       g2hjWW,g2hjZZ,g2hjZga,g2hjgaga,g2hjgg,g2hjggZ,
     &       g2hjhiZ,BR_hjhihi,BR_hjinvisible

      nHzero=1
      nHplus=0
         
      call initialize_HiggsBounds(nHzero,nHplus,'LandH')
c         call initialize_HiggsBounds(nHzero,nHplus,'onlyL')
c         call initialize_HiggsBounds(nHzero,nHplus,'onlyH')
c         call initialize_HiggsBounds(nHzero,nHplus,'onlyP')

      open(21,file='example-SM-results.dat')
      open(20,file='example-4thGen-results.dat')

      do i=90,600,1
        Mh=dble(i)
        SMGammaTotal=SMGamma_h(Mh)

! SMGamma_h(Mh), SMBR_Hgg(Mh), SMBR_Hgg(Mh) are set to -1 if called
! with Mh out of range [0.8 GeV, 500 GeV]. The calculation is then bypassed.

        if(.not. (SMGammaTotal .lt. 0)) then
! SM 
          GammaTotal=SMGammaTotal
c CP even
          g2hjss_s=1d0
          g2hjss_p=0d0
          g2hjcc_s=1d0
          g2hjcc_p=0d0
          g2hjbb_s=1d0
          g2hjbb_p=0d0
          g2hjtt_s=1d0
          g2hjtt_p=0d0         
          g2hjmumu_s=1d0
          g2hjmumu_p=0d0  
          g2hjtautau_s=1d0
          g2hjtautau_p=0d0

          g2hjWW=1d0
          g2hjZZ=1d0
          g2hjZga=1d0
          g2hjgaga=1d0
          g2hjgg=1d0
          g2hjggZ=1d0
          g2hjhiZ=0d0
          BR_hjhihi=0d0
          BR_hjinvisible=0d0

          call HiggsBounds_neutral_input_effC(Mh,GammaTotal,
     &    g2hjss_s,g2hjss_p,g2hjcc_s,g2hjcc_p,g2hjbb_s,g2hjbb_p,
     &    g2hjtt_s,g2hjtt_p,
     &    g2hjmumu_s,g2hjmumu_p,g2hjtautau_s,g2hjtautau_p,
     &    g2hjWW,g2hjZZ,g2hjZga,g2hjgaga,g2hjgg,g2hjggZ,
     &    g2hjhiZ, BR_hjinvisible,BR_hjhihi)

          call run_HiggsBounds(HBresult,chan,obsratio,ncombined)
          write(21,*) Mh,HBresult,chan,obsratio

! 4th Generation Model (using approximate enhancement factor 9 for
!                       the hgg coupling)
          SMGamma_hgg=SMBR_Hgg(Mh)*SMGammaTotal
          Gamma_hgg=9d0*SMGamma_hgg
          FGMGammaTotal=SMGammaTotal-SMGamma_hgg+Gamma_hgg
          GammaTotal=FGMGammaTotal

          g2hjgaga=1d0/9d0*FGMGammaTotal/SMGammaTotal
          g2hjgg=9d0

          call HiggsBounds_neutral_input_effC(Mh,GammaTotal,
     &    g2hjss_s,g2hjss_p,g2hjcc_s,g2hjcc_p,g2hjbb_s,g2hjbb_p,
     &    g2hjtt_s,g2hjtt_p,
     &    g2hjmumu_s,g2hjmumu_p,g2hjtautau_s,g2hjtautau_p,
     &    g2hjWW,g2hjZZ,g2hjZga,g2hjgaga,g2hjgg,g2hjggZ,
     &    g2hjhiZ, BR_hjinvisible,BR_hjhihi)

          call run_HiggsBounds(HBresult,chan,obsratio,ncombined)
          write(20,*) Mh,HBresult,chan,obsratio
        endif
      enddo

      close(21)
      close(20)

      call finish_HiggsBounds

      end
