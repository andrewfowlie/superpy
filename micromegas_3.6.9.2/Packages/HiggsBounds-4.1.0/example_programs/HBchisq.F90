!--------------------------------------------------------------------------------------
! This example program is part of HiggsBounds (TS 07/04/2013).
!--------------------------------------------------------------------------------------
program HBchisq
! This example program scans the SM Higgs mass from 100 GeV to 120 GeV and evaluates
! the usual HiggsBounds 95% C.L. limit as well as the chi^2 value from the LEP exclusions.
!--------------------------------------------------------------------------------------
 use theory_colliderSfunctions
 implicit none

  integer :: nHzero, nHplus
  integer :: HBresult,chan,ncombined,chan2
  integer,parameter :: fileid=78, fileid2=80
  double precision :: obsratio,theory_uncertainty_1s,chisq_withouttheory,chisq_withtheory
  double precision :: SMGamma_h,SMGammaTotal
  double precision :: Mh,g2hjss_s,g2hjss_p,g2hjcc_s,g2hjcc_p, &
&                    g2hjbb_s,g2hjbb_p,g2hjtt_s,g2hjtt_p, &
&                    g2hjmumu_s,g2hjmumu_p,g2hjtautau_s,g2hjtautau_p, &
&                    g2hjWW,g2hjZZ,g2hjZga,g2hjgaga,g2hjgg,g2hjggZ,	&
&                    g2hjhiZ,BR_hjhihi,BR_hjinvisible
  character(len=100)::filename
  integer :: i
  
  nHzero=1
  nHplus=0
  theory_uncertainty_1s=1.5D0

!--Setting up the output
 filename='HBchisq-output.dat'
 open(fileid,file=filename)
 write(fileid,*) '# mh   HBres   channel   obsratio   ncombined   chisq_wo_uncertainty ',&
& '  chisq_w_uncertainty channel(chisq)' 
 write(fileid,*) '#--------------------------------------------------------------------',&
& '----------------------#'

!---- Initialize HiggsBounds with the LEP chisq tables ----!
 call initialize_HiggsBounds_chisqtables
!---- Initialize HiggsBounds with only LEP results ----!
 call initialize_HiggsBounds(nHzero, nHplus, "onlyL")

 do i=1,101
   Mh = 100.0D0 + (i-1)*0.2D0

   SMGammaTotal=SMGamma_h(Mh)

! SMGamma_h(Mh), SMBR_Hgg(Mh), SMBR_Hgg(Mh) are set to -1 if called
! with Mh out of range [0.8 GeV, 500 GeV]. The calculation is then bypassed.
   if(.not. (SMGammaTotal .lt. 0)) then
    g2hjss_s=1.0d0
    g2hjss_p=0.0d0
    g2hjcc_s=1.0d0
    g2hjcc_p=0.0d0
    g2hjbb_s=1.0d0
    g2hjbb_p=0.0d0
    g2hjtt_s=1.0d0
    g2hjtt_p=0.0d0         
    g2hjmumu_s=1.0d0
    g2hjmumu_p=0.0d0  
    g2hjtautau_s=1.0d0
    g2hjtautau_p=0.0d0
    g2hjWW=1.0d0
    g2hjZZ=1.0d0
    g2hjZga=1d0
    g2hjgg=1.0d0
    g2hjggZ=1d0
    g2hjhiZ=0d0
    g2hjgaga=1.0d0
    BR_hjhihi=0d0
    BR_hjinvisible=0d0
      
    call HiggsBounds_neutral_input_effC(Mh,SMGammaTotal, &
     &    g2hjss_s,g2hjss_p,g2hjcc_s,g2hjcc_p,g2hjbb_s,g2hjbb_p, &
     &    g2hjtt_s,g2hjtt_p, &
     &    g2hjmumu_s,g2hjmumu_p,g2hjtautau_s,g2hjtautau_p, &
     &    g2hjWW,g2hjZZ,g2hjZga,g2hjgaga,g2hjgg,g2hjggZ, &
     &    g2hjhiZ, BR_hjinvisible,BR_hjhihi)

    call run_HiggsBounds_classic( HBresult, chan, obsratio, ncombined )
    
    call HB_calc_stats(theory_uncertainty_1s,chisq_withouttheory,chisq_withtheory,chan2)

    write(fileid,*) Mh, HBresult, chan, obsratio, ncombined, &
&                   chisq_withouttheory,chisq_withtheory, chan2 

   endif

 enddo
 
 close(fileid)

  call finish_HiggsBounds_chisqtables
  call finish_HiggsBounds

end program HBchisq