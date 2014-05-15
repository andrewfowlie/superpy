!--------------------------------------------------------------------------------------
program HS_scale_uncertainties
! This example shows how the user can assign scale factors to the peak observables,
! which scale the rate (mu) uncertainty. These factors can be assigned individually
! to the observable specified by its ID. Here, we use a universal scale factor to keep
! things simple, and scan over its values.
!   We first scale only the experimental rate uncertainties (which can be motivated by
! having a larger data sample, for instance), while leaving the theoretical rate 
! uncertainties fixed. In the second step, we scale only the theory uncertainties.
! Finally, we scale both experimental and theoretical rate uncertainties by the same
! factor. The results is stored in textfiles which can be plotted.
! As model we use the SM with a fixed Higgs mass. We use the effective couplings input.
!
! This example program is part of HiggsSignals (TS 28/01/2013).
!--------------------------------------------------------------------------------------
 use theory_colliderSfunctions
 use usefulbits, only : vsmall
 use io, only : get_number_of_observables, get_ID_of_peakobservable
 implicit none

 integer :: nHzero, nHplus, ndf, i, j, k, ii, jj
 double precision :: obsratio, mass, Pvalue, Chisq, mu, Chisq_mu, Chisq_mh
 double precision :: dCS(5),dBR(5)
 double precision :: SMGammaTotal, SMGamma_h
 double precision :: Mh,GammaTotal,g2hjss_s,g2hjss_p,g2hjcc_s,g2hjcc_p,					&
     &       g2hjbb_s,g2hjbb_p,g2hjtt_s,g2hjtt_p,										&
     &       g2hjmumu_s,g2hjmumu_p,g2hjtautau_s,g2hjtautau_p,							&
     &       g2hjWW,g2hjZZ,g2hjZga,g2hjgaga,g2hjgg,g2hjggZ,								&
     &       g2hjhiZ,BR_hjhihi,BR_hjinvisible
 character(len=100)::filename
 double precision :: dm
 integer		  :: pdf
 double precision :: sf
 integer :: ntotal, npeakmu, npeakmh, nmpred, nanalyses, ID

 
 nHzero=1
 nHplus=0

!--Set the (relative!) rate uncertainties for your model:
!  dCS(1) - singleH				dBR(1) - gamma gamma
!  dCS(2) - VBF					dBR(2) - W W
!  dCS(3) - HW					dBR(3) - Z Z
!  dCS(4) - HZ					dBR(4) - tau tau
!  dCS(5) - ttH					dBR(5) - b bbar
!--Here, we first set them to the SM. Taken SM values from
!--https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CrossSections
 dCS = (/ 0.147D0, 0.028D0, 0.037D0, 0.051D0, 0.12D0 /)
 dBR = (/ 0.054D0, 0.048D0, 0.048D0, 0.061D0, 0.028D0 /)
 
!--Enter the Higgs mass and its theory uncertainty here: 
 Mh = 125.8D0
 pdf = 2
!---- Initialize HiggsSignals and pass the name of the experimental analysis folder  ----!
 call initialize_HiggsSignals(nHzero,nHplus,"latestresults")
!---- Set the output level (0: silent, 1: screen output, etc...)				 	 ----!
 call setup_output_level(1)
!---- Set the Higgs mass parametrization (1: box, 2:Gaussian, 3:Gox+gaussian)	 	 ----!
 call setup_pdf(pdf) 
 !---- Set the assignment range for the peak-centered method (optional)				 ----! 
 call setup_assignmentrange(2.0D0)
 !---- Disable the use of CS and BR uncertainties via provided covariance matrices	 ----!
 !     This is needed in this simple example (since HiggsSignals-1.1.0) because we
 !     scale the (naively) estimated maximal rate uncertainties. A more sophisticated
 !     treatment of scaled rate uncertainties requires changing the covariance matrix
 !     for the XS and BRs for each projected point.
 call setup_correlated_rate_uncertainties(0)
  
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
  g2hjmumu_s=1d0
  g2hjmumu_p=0.0d0  
  g2hjtautau_s=1d0
  g2hjtautau_p=0d0
  g2hjWW=1.0d0
  g2hjZZ=1.0d0
  g2hjZga=1.0d0
  g2hjgg=1.0d0
  g2hjggZ=1.0d0
  g2hjhiZ=0d0
  g2hjgaga=1d0
  BR_hjhihi=0d0
  BR_hjinvisible=0d0      

!-Here, we set the theoretical rate uncertainties to their original SM values.
  call setup_rate_uncertainties(dCS, dBR)

!-We first perform a dummy run to obtain the number of observables and their IDs.
  call HiggsBounds_neutral_input_effC(Mh,SMGammaTotal,									&
     &    g2hjss_s,g2hjss_p,g2hjcc_s,g2hjcc_p,g2hjbb_s,g2hjbb_p,						&
     &    g2hjtt_s,g2hjtt_p,															&
     &    g2hjmumu_s,g2hjmumu_p,g2hjtautau_s,g2hjtautau_p,								&
     &    g2hjWW,g2hjZZ,g2hjZga,g2hjgaga,g2hjgg,g2hjggZ,								&
     &    g2hjhiZ, BR_hjinvisible,BR_hjhihi)
   call run_HiggsSignals(1, Chisq_mu, Chisq_mh, Chisq, ndf, Pvalue)
   
!-Get the number of the peak-observables (Don't care about ntotal,npeakmh,nmpred,nanalyses)
  call get_number_of_observables(ntotal, npeakmu, npeakmh, nmpred, nanalyses)
 
!---- For the scan, be silent.													 	 ----!
 call setup_output_level(0)
  
!-1) Scale only the experimental mu uncertainties and leave the theory uncertainties
!    as original.
   write(*,*) "Scaling only exp. rate uncertainties..."

!--Setting up the output
  filename='results/scaling_dmu_exp.dat'
  open(21,file=filename)  
  write(21,'(A10,1X,A10)') "# Scale-F."," chi^2"

  do i=1,100
!--In this example, the individual rate uncertainties are scaled by the same factor:  
  sf = 1.0 - (i-1)*0.01

!-We now loop over all peak observables (known from the dummy run) and obtain for each
!-observable its ID. In the same step, we assign the scale factor sf to this observable.
  do ii=1,npeakmu
   call get_ID_of_peakobservable(ii, ID)
   call assign_rate_uncertainty_scalefactor_to_peak(ID, sf)
  enddo

!!-Set the HiggsSignals input (again!)
   call HiggsBounds_neutral_input_effC(Mh,SMGammaTotal,							&
     &    g2hjss_s,g2hjss_p,g2hjcc_s,g2hjcc_p,g2hjbb_s,g2hjbb_p,						&
     &    g2hjtt_s,g2hjtt_p,															&
     &    g2hjmumu_s,g2hjmumu_p,g2hjtautau_s,g2hjtautau_p,								&
     &    g2hjWW,g2hjZZ,g2hjZga,g2hjgaga,g2hjgg,g2hjggZ,								&
     &    g2hjhiZ, BR_hjinvisible,BR_hjhihi)
!-Run HiggsSignals with scaled uncertainties for mu. Note that after the run, the
!-uncertainties are set back to their original values.
   call run_HiggsSignals(1, Chisq_mu, Chisq_mh, Chisq, ndf, Pvalue)
! This will collect the main HiggsSignals results together into one file
    write(21,'(2F20.4)')sf,Chisq
  enddo
  close(21)

!-2) Scale only the theoretical mu uncertainties and leave the experimental uncertainties
!    as original. We scale only the uncertainties of the production modes.
   write(*,*) "Scaling only th. rate uncertainties..."

!--Setting up the output
  filename='results/scaling_dmu_th.dat'
  open(21,file=filename)  
  write(21,'(A10,1X,A10)') "# Scale-F."," chi^2"

  do i=1,100
!--In this example, the individual rate uncertainties are scaled by the same factor:  
  sf = 1.0 - (i-1)*0.01

!-Here, we set the theoretical rate uncertainties to their original SM values.
  call setup_rate_uncertainties(dCS*sf, dBR*sf)

!!-Set the HiggsSignals input (again!)
   call HiggsBounds_neutral_input_effC(Mh,SMGammaTotal,							&
     &    g2hjss_s,g2hjss_p,g2hjcc_s,g2hjcc_p,g2hjbb_s,g2hjbb_p,						&
     &    g2hjtt_s,g2hjtt_p,															&
     &    g2hjmumu_s,g2hjmumu_p,g2hjtautau_s,g2hjtautau_p,								&
     &    g2hjWW,g2hjZZ,g2hjZga,g2hjgaga,g2hjgg,g2hjggZ,								&
     &    g2hjhiZ, BR_hjinvisible,BR_hjhihi)
!-Run HiggsSignals with scaled uncertainties for mu. Note that after the run, the
!-uncertainties are set back to their original values.
   call run_HiggsSignals(1, Chisq_mu, Chisq_mh, Chisq, ndf, Pvalue)
! This will collect the main HiggsSignals results together into one file
    write(21,'(2F20.4)')sf,Chisq
  enddo
  close(21)

!-3) Scale both experimental and theoretical mu uncertainties by same scale factor.
   write(*,*) "Scaling both exp. and th. rate uncertainties..."

!--Setting up the output
  filename='results/scaling_dmu_both.dat'
  open(21,file=filename)  
  write(21,'(A10,1X,A10)') "# Scale-F."," chi^2"

  do i=1,100 
!--In this example, the individual rate uncertainties are scaled by the same factor:  
  sf = 1.0 - (i-1)*0.01

!-We now loop over all peak observables (known from the dummy run) and obtain for each
!-observable its ID. In the same step, we assign the scale factor sf to this observable.
  do ii=1,npeakmu
   call get_ID_of_peakobservable(ii, ID)
   call assign_rate_uncertainty_scalefactor_to_peak(ID, sf)
  enddo

!-Here, we set the theoretical rate uncertainties to their original SM values.
  call setup_rate_uncertainties(dCS*sf, dBR*sf)

!!-Set the HiggsSignals input (again!)
   call HiggsBounds_neutral_input_effC(Mh,SMGammaTotal,							&
     &    g2hjss_s,g2hjss_p,g2hjcc_s,g2hjcc_p,g2hjbb_s,g2hjbb_p,						&
     &    g2hjtt_s,g2hjtt_p,															&
     &    g2hjmumu_s,g2hjmumu_p,g2hjtautau_s,g2hjtautau_p,								&
     &    g2hjWW,g2hjZZ,g2hjZga,g2hjgaga,g2hjgg,g2hjggZ,								&
     &    g2hjhiZ, BR_hjinvisible,BR_hjhihi)
!-Run HiggsSignals with scaled uncertainties for mu. Note that after the run, the
!-uncertainties are set back to their original values.
   call run_HiggsSignals(1, Chisq_mu, Chisq_mh, Chisq, ndf, Pvalue)
   
! This will collect the main HiggsSignals results together into one file
    write(21,'(2F20.4)')sf,Chisq
  enddo
  close(21)


 endif
  
 call finish_HiggsSignals

end program HS_scale_uncertainties