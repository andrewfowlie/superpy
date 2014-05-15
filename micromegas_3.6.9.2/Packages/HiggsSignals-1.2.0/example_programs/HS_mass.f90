!--------------------------------------------------------------------------------------
! This example program is part of HiggsSignals (TS 26/09/2013).
!--------------------------------------------------------------------------------------
program HS_mass
! In this example the peak-centered chi^2 method is applied to a SM-like Higgs boson
! (overall signal strength scale factor mu) within the mass range 110 - 140 GeV. 
! All three mass pdf choices are considered. Theoretical mass uncertainties and
! assignment range can be changed.
!--------------------------------------------------------------------------------------
 use theory_colliderSfunctions
 use usefulbits, only : vsmall
 use pc_chisq, only : print_cov_mh_to_file,print_cov_mu_to_file,print_inverse_cov_mh_to_file,&
&                    get_peakchi2, print_corr_mu_to_file 
 use io, only : get_number_of_observables,get_ID_of_peakobservable,get_peakinfo_from_HSresults
 implicit none

 integer :: nHzero, nHplus, ndf, i, j, k, ii, jj
 double precision :: obsratio, mass, Pvalue, Chisq, mu, Chisq_mu, Chisq_mh, Lambda
 double precision :: dCS(5),dBR(5)
 double precision :: SMGammaTotal
 double precision :: scale_bbh, scale_ggh, dggh, dbbh
 double precision :: mh,GammaTotal,g2hjss_s,g2hjss_p,g2hjcc_s,g2hjcc_p, &
&                    g2hjbb_s,g2hjbb_p,g2hjtt_s,g2hjtt_p, &
&                    g2hjmumu_s,g2hjmumu_p,g2hjtautau_s,g2hjtautau_p, &
&                    g2hjWW,g2hjZZ,g2hjZga,g2hjgaga,g2hjgg,g2hjggZ,	&
&                    g2hjhiZ,BR_hjhihi,BR_hjinvisible
 character(len=100)::filename
 double precision :: dm
 integer		  :: pdf
 integer :: ntotal, npeakmu, npeakmh, nmpred, nanalyses, ID, domH, nHcomb, Nassigned
 double precision :: mupred
 double precision, allocatable :: csqmu(:),csqmh(:),csqmax(:),csqtot(:)
 integer, allocatable :: ncomb(:)
!-HiggsBounds internal functions to obtain SM branching ratios 
 double precision :: SMBR_Htoptop,SMBR_Hss, SMBR_Hcc, SMBR_Hbb, SMBR_Hmumu, SMBR_Htautau,&
 &                   SMBR_HWW, SMBR_HZZ, SMBR_HZgam, SMBR_Hgamgam, SMBR_Hgg,SMGamma_h

 nHzero=1
 nHplus=0

!---- Initialize HiggsSignals and pass the name of the experimental analysis folder ----!
 call initialize_HiggsSignals(nHzero,nHplus,"latestresults")
!---- Set the output level (0: silent, 1: screen output, 2: even more output,...) ----!
 call setup_output_level(0)

!--Enter the Higgs mass and its theory uncertainty here: 
 dm = 2.0D0
 mu = 1.0D0
 Lambda = 1.0D0
 
!---- Pass the Higgs mass uncertainty to HiggsSignals ----!
 call HiggsSignals_neutral_input_MassUncertainty(dm)
!---- Set the assignment range for the peak-centered method (optional)				 ----! 
!     This can be done either to all observables or only to the 
!     mass-sensitive observables, which contribute to the Higgs mass chi^2
! call setup_assignmentrange(Lambda)
 call setup_assignmentrange_massobservables(Lambda)

 call setup_correlations(1,1) ! mu, mass

 do pdf=1,3
!---- Set the Higgs mass parametrization (1: box, 2:gaussian, 3:box+gaussian) ----!
  call setup_pdf(pdf) 
  select case(pdf)
   case(1)
    filename='results/HS_mass_pdf1.dat'
   case(2)
    filename='results/HS_mass_pdf2.dat'
   case(3)
    filename='results/HS_mass_pdf3.dat'
   case default
  end select
  
  open(21,file=filename)
  write(21,*) '# mh   dmh     Chisq_mu    Chisq_mh    Chisq    ndf' 
  write(21,*) '#----------------------------------------------------'

  do j=1,301 !181,181!
   mh = 110.0D0 +(j-1)*0.1D0

   SMGammaTotal=SMGamma_h(Mh)

! SMGamma_h(Mh), SMBR_Hgg(Mh), SMBR_Hgg(Mh) are set to -1 if called
! with Mh out of range [0.8 GeV, 500 GeV]. The calculation is then bypassed.
   if(.not. (SMGammaTotal .lt. 0)) then
    g2hjss_s=1.0d0*mu
    g2hjss_p=0.0d0*mu
    g2hjcc_s=1.0d0*mu
    g2hjcc_p=0.0d0*mu
    g2hjbb_s=1.0d0*mu
    g2hjbb_p=0.0d0*mu
    g2hjtt_s=1.0d0*mu
    g2hjtt_p=0.0d0*mu         
    g2hjmumu_s=1.0d0*mu
    g2hjmumu_p=0.0d0*mu
    g2hjtautau_s=1.0d0*mu
    g2hjtautau_p=0.0d0*mu
    g2hjWW=1.0d0*mu
    g2hjZZ=1.0d0*mu
    g2hjZga=1d0*mu
    g2hjgg=1.0d0*mu
    g2hjggZ=1d0*mu
    g2hjhiZ=0d0*mu
    g2hjgaga=1.0d0*mu
    BR_hjhihi=0d0
    BR_hjinvisible=0d0
      
!----Calculate the new total decay width:
    GammaTotal = SMGammaTotal*(1 + &
	&  	            (g2hjWW - 1)*SMBR_HWW(Mh)+(g2hjZZ - 1)*SMBR_HZZ(Mh) + &
	&               (g2hjgg - 1)*SMBR_Hgg(Mh)+(g2hjtt_s - 1)*SMBR_Htoptop(Mh)+ &
	&               (g2hjbb_s - 1)*SMBR_Hbb(Mh)+(g2hjtautau_s - 1)*SMBR_Htautau(Mh)+ &
	&               (g2hjss_s - 1)*SMBR_Hss(Mh)+(g2hjcc_s - 1)*SMBR_Hcc(Mh)+ &
	&               (g2hjZga - 1)*SMBR_HZgam(Mh)+(g2hjmumu_s - 1)*SMBR_Hmumu(Mh)+ &		
	&               (g2hjgaga - 1)*SMBR_Hgamgam(Mh)	)


    call HiggsBounds_neutral_input_effC(mh,GammaTotal, &
     &    g2hjss_s,g2hjss_p,g2hjcc_s,g2hjcc_p,g2hjbb_s,g2hjbb_p, &
     &    g2hjtt_s,g2hjtt_p, &
     &    g2hjmumu_s,g2hjmumu_p,g2hjtautau_s,g2hjtautau_p, &
     &    g2hjWW,g2hjZZ,g2hjZga,g2hjgaga,g2hjgg,g2hjggZ, &
     &    g2hjhiZ, BR_hjinvisible,BR_hjhihi)

    call run_HiggsSignals( 1, Chisq_mu, Chisq_mh, Chisq, ndf, Pvalue)

    call get_number_of_observables(ntotal, npeakmu, npeakmh, nmpred, nanalyses)

    allocate(csqmu(npeakmu),csqmh(npeakmu),csqmax(npeakmu),csqtot(npeakmu),ncomb(npeakmu))

    Nassigned=0
    do ii=1,npeakmu
     call get_ID_of_peakobservable(ii, ID)
     call get_peakinfo_from_HSresults(ID, mupred, domH, nHcomb)
     ncomb(ii)=nHcomb
     call get_peakchi2(ID, csqmu(ii), csqmh(ii), csqmax(ii), csqtot(ii))
     Nassigned=Nassigned+nHcomb
    enddo 

    deallocate(csqmu,csqmh,csqmax,csqtot,ncomb)

    write(21,*) mh,dm,mu,Chisq_mu,Chisq_mh,Chisq,Nassigned,ndf,Lambda

   endif
  enddo
  close(21)
 enddo
 

 write(*,*) "Finishing HiggsSignals..."
 call finish_HiggsSignals

 end program HS_mass
