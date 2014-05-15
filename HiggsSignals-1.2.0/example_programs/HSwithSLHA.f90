!--------------------------------------------------------------------------------------
program HSwithSLHA
!--------------------------------------------------------------------------------------
! This example program is part of HiggsSignals (TS 29/01/2013).
!
! In this example we demonstrate how HiggsSignals can be run on SLHA files. The SLHA
! file has to contain the two HiggsBounds SLHA input blocks
!	HiggsBoundsInputHiggsCouplingsBosons
!	HiggsBoundsInputHiggsCouplingsFermions
! (see HiggsBounds (version 3 or more) manual for more details)
!
! Run with
!   ./HSwithSLHA npoints <stem>
! where npoints is the number of parameter points you would like to
! look at and each parameter point has a corresponding SLHA file
! e.g. the corresponding SLHA for the 5th point should be found at
!   <stem>.5
!
! Output:
! The HiggsSignals SLHA output blocks will be added to each SLHA file.
! The results are summarized in an additional file <stem>-fromHS
!
! We furthermore demonstrate how to get the signal rates directly from HiggsSignals
! after a successful run.
!--------------------------------------------------------------------------------------
  use io, only : HiggsSignals_SLHA_output, get_peakinfo_from_HSresults
  use pc_chisq, only : print_cov_mu_to_file, print_peaks_to_file
  implicit none
  integer :: nH,nHplus,ndf
  double precision :: Chisq, Chisq_mu, Chisq_mh, Pvalue
  double precision :: dCS(5),dBR(5)
  double precision :: R_H_WW, R_H_ZZ, R_H_gaga, R_H_tautau, R_H_bb, R_VH_bb, 			&
&                     totalrate, ggf_rate, hgaga_rate, gghgg_rate
  double precision, allocatable :: dMh(:), masses(:), dmasses(:)
  integer :: i,npoints
  integer,parameter :: fileid=78, fileid2=79
  character(len=8) :: istring
  character(len=300) :: inputfilename,outputfilename
  character(len=300) :: stem
  character(LEN=300) :: temp, tmpstring
  integer :: number_args, ios
!--------------------------------------------------------------------------------------
  nH=3
  nHplus=1

  allocate(dMh(nH),masses(nH),dmasses(nH))
!--Give estimates on (relative!) systematic uncertainties for:
!  dCS(1) - singleH				dBR(1) - gamma gamma
!  dCS(2) - VBF					dBR(2) - W W
!  dCS(3) - HW					dBR(3) - Z Z
!  dCS(4) - HZ					dBR(4) - tau tau
!  dCS(5) - ttH					dBR(5) - b bbar
!--Enter now the rate uncertainties for your model (typical MSSM estimates):
  dCS = (/ 0.20D0, 0.028D0, 0.037D0, 0.051D0, 0.12D0 /)
  dBR = (/ 0.054D0, 0.048D0, 0.048D0, 0.061D0, 0.028D0 /)
!--n.b. have to set theoretical uncertainties on Higgs masses dMh (in GeV) for h,H,A:
! These are default values. In SLHA mode, these values are read in from the Block "DMASS".
  dMh = (/ 2.0D0, 2.0D0, 0.0D0 /)  
!----------------------------------- preprocessing --------------------------------------!
  number_args = IARGC() 
  if( number_args .ne. 2)then
   stop "Incorrect number of arguments given to HSwithSLHA"
  endif
  ! Read arguments into text strings.
  i=1
  temp=""
  call GETARG(i,temp)
  read(temp,*) npoints
  i=i+1  
  temp=""
  call GETARG(i,temp)
  stem = ""
  stem = trim(temp)
!------------------------------ HiggsSignals options ------------------------------------!
!---- Initialize HiggsSignals and pass the name of the experimental analysis folder  ----! 
 call initialize_HiggsSignals(nH,nHplus,"latestresults")  
!---- Set the output level (0: silent, 1: screen output, 2: even more output,...) 	 ----!
 call setup_output_level(1)
!---- Set the Higgs mass parametrization (1: box, 2:gaussian, 3:box+gaussian)	 	 ----!
 call setup_pdf(2) 
!---- Set the assignment range for the peak-centered method (optional)				 ----! 
 call setup_assignmentrange_massobservables(2.0D0)
!---- If the mass-centered chi^2 method is used, can specify the dm_theory treatment ----!
!     1: mass variation, 2: smearing of mu-plot with mass pdf  
! call setup_mcmethod_dm_theory(1)
!---- Use symmetric rate errors? (0: original(default), 1: averaged-symmetrical)     ----!
! call setup_symmetricerrors(0) 
!---- Are SM rate uncertainties included in signal strength measurements?            ----!
!    (0: no, 1: yes(default) )
! call setup_thu_observables(1)
!---- Allow anti-correlated signal strength measurements? (0: no, 1: yes(default) )  ----!
! call setup_anticorrelations_in_mu(1)
!----


  outputfilename=trim(adjustl(stem))//'-fromHS'
  open(fileid, file=trim(outputfilename))

  do i=1,npoints
   if(i.gt.99999999)stop'need to increase the size of istring in HSwithSLHA'
   write(istring,'(I8)')i
   inputfilename=trim(adjustl(stem))//'.'//trim(adjustl(istring))
!--Test if input file exists and is non-empty
   open(fileid2, file=inputfilename, form='formatted')
   read(fileid2,'(A)',iostat=ios) tmpstring

   if(ios.eq.0) then
    close(fileid2)  
    
!-------------------------------- HiggsSignals run --------------------------------------!    
!---- Feed HiggsSignals with the the model input using HiggsBounds subroutine	 	 ----!
    call HiggsBounds_input_SLHA(inputfilename)
!---- Checking the Higgs mass uncertainty           							 	 ----!
! The theoretical Higgs mass uncertainties are read in from the SLHA Block "DMASS"
! in the call of the subroutine HiggsBounds_input_SLHA. If the block "DMASS" is absent,
! they are set to zero. If the user wants to change the values obtained from the SLHA
! file, he/she can call the subroutine
! 	HiggsSignals_neutral_input_MassUncertainty
! AFTER reading in the SLHA file.
! Here, we only print out the values which have been stored already in 
! HiggsBounds/HiggsSignals. If their sum is <= zero, the Block DMASS was probably absent
! and we set the default values as specified above.
!
    call get_neutral_Higgs_masses(masses, dmasses)
    write(*,*) "Neutral Higgs boson mass spectrum (from SLHA): "
    write(*,'(2X,A,F6.2,A,F4.2)') "mass(h0) = ",masses(1)," +- ",dmasses(1)
    write(*,'(2X,A,F6.2,A,F4.2)') "mass(H0) = ",masses(2)," +- ",dmasses(2)
    write(*,'(2X,A,F6.2,A,F4.2)') "mass(A0) = ",masses(3)," +- ",dmasses(3)
    if(sum(dmasses).le.0.0D0) then
     write(*,*) "BLOCK DMASS not found, changing mass uncertainies to ", dMh
 	 call HiggsSignals_neutral_input_MassUncertainty(dMh)
 	endif 
!---- Set the production and decay rate uncertainties for the model			 	 	 ----!
    call setup_rate_uncertainties(dCS, dBR)
!---- Run HiggsSignals (runmode: 1 (peak-centered), 2 (mass-centered), 3 (both)
    call run_HiggsSignals(1, Chisq_mu, Chisq_mh, Chisq, ndf, Pvalue)
!----------------------------- HiggsSignals output --------------------------------------!    
!---- Attach HiggsSignals SLHA blocks to SLHA file									 ----!
!   integer argument gives level of details:
!   0 : writes only HiggsSignalsResults block
!   else : writes all blocks
    call HiggsSignals_SLHA_output(1)
!---- Now, some examples of how to read out the signal rates. Note, that these subroutines
!     have to be called after run_HiggsSignals
!
!---- Get signal-rate ratios (without efficiencies) for lightest Higgs boson and LHC8----!
    call get_Rvalues(1, 3, R_H_WW, R_H_ZZ, R_H_gaga, R_H_tautau, R_H_bb, R_VH_bb)
!---- Get the total signal rate (without efficiencies)					 			 ----!
	call get_rates(1,3,25,(/11,12,13,14,15,21,22,23,24,25,31,32,33,34,35,				&
&							41,42,43,44,45,51,52,53,54,55/),totalrate)
!---- Get the gluon gluon fusion rate (without efficiencies)			 			 ----!
	call get_rates(1,3,1,(/10/),ggf_rate)
!---- Get the H -> gamma gamma rate (without efficiencies)				 			 ----!
	call get_rates(1,3,1,(/01/),hgaga_rate)
!	NEW SINCE HiggsSignals-1.1.0: more decay modes accessible via get_rates:
!	  Decay mode ID (Final state):  6 (Zgamma), 7 (cc), 8 (mumu), 9 (gg)
!---- Get the gg->H->gg (without efficiencies)				 			             ----!
	call get_rates(1,3,1,(/19/),gghgg_rate)

    write(*,'(A,F10.4)') "R_H_WW     = ", R_H_WW
    write(*,'(A,F10.4)') "R_H_ZZ     = ", R_H_ZZ
    write(*,'(A,F10.4)') "R_H_gaga   = ", R_H_gaga
    write(*,'(A,F10.4)') "R_H_tautau = ", R_H_tautau
    write(*,'(A,F10.4)') "R_H_bb     = ", R_H_bb
    write(*,'(A,F10.4)') "R_VH_bb    = ", R_VH_bb
    write(*,'(A,F10.4)') "totalrate  = ", totalrate
    write(*,'(A,F10.4)') "ggf_rate   = ", ggf_rate
    write(*,'(A,F10.4)') "h->gaga    = ", hgaga_rate
    write(*,'(A,F10.4)') "gg->h->gg  = ", gghgg_rate

!--This will collect the main HiggsSignals results together into one file
    write(fileid,*)i,Pvalue,Chisq,Chisq_mu,Chisq_mh,ndf
   else
    close(fileid2)     
    call system("rm -f "//inputfilename)    
   endif 

  enddo
  
  close(fileid)
! call print_peaks_to_file
! call print_cov_mu_to_file

 call finish_HiggsSignals
    
end program HSwithSLHA