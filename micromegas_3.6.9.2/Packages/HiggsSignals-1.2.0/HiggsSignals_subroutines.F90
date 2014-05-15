!------------------------------------------------------------
! This file is part of HiggsSignals (TS 03/03/2013).
!------------------------------------------------------------
subroutine initialize_HiggsSignals_latestresults(nHiggsneut,nHiggsplus)
!------------------------------------------------------------
! Wrapper subroutine to intitialize HiggsSignals with the experimental
! dataset "latestresults", avoiding to specify this via a string argument.
!------------------------------------------------------------
 implicit none
 !--------------------------------------input
 integer,intent(in) :: nHiggsneut
 integer,intent(in) :: nHiggsplus
 character(LEN=13) :: Expt_string
 
 Expt_string = "latestresults"

 call initialize_HiggsSignals(nHiggsneut,nHiggsplus,Expt_string)
 
end subroutine initialize_HiggsSignals_latestresults
!------------------------------------------------------------
subroutine initialize_HiggsSignals(nHiggsneut,nHiggsplus,Expt_string)
!------------------------------------------------------------
! This the first HiggsSignals subroutine that should be called
! by the user.
! It calls subroutines to read in the tables of Standard Model 
! decay and production rates from HiggsBounds, sets up the 
! experimental data from Tevatron and LHC, allocate arrays, etc.
! Arguments (input):
!   * nHiggs = number of neutral Higgs in the model 
!   * nHiggsplus = number of singly, positively charged Higgs in the model
!   * Expt_string = name of experimental dataset to be used
!------------------------------------------------------------
 use usefulbits, only : np,Hneut,Hplus,Chineut,Chiplus,debug,inputmethod,&
  &   inputsub,theo,whichanalyses,just_after_run,&
  &   file_id_debug1,file_id_debug2,allocate_if_stats_required
 use usefulbits_HS, only : HiggsSignals_info, nanalys, eps, Exptdir, obs
 use datatables, only: setup_observables
 use input, only : check_number_of_particles,check_whichanalyses
 use io, only : setup_input_for_hs, setup_output_for_hs
 use theory_BRfunctions, only : setup_BRSM, BRSM

#if defined(NAGf90Fortran)
 use F90_UNIX_IO, only : flush
#endif

 implicit none
 !--------------------------------------input
 integer,intent(in) :: nHiggsneut
 integer,intent(in) :: nHiggsplus
 character(LEN=*), intent(in) :: Expt_string
 !-----------------------------------internal
 integer :: i
 !----------------------------------parameter
 eps=5.0D0
 np(Hneut)=nHiggsneut

 np(Hplus)=nHiggsplus
 
 Exptdir = Expt_string
 
 np(Chineut)=0! not considering bounds on neutralinos here
 np(Chiplus)=0! not considering bounds on charginos here
        
 debug=.False.

 select case(whichanalyses)
  case('onlyL')
   whichanalyses='LandH'
  case('onlyH','onlyP','list ','LandH')
  case default
   whichanalyses='onlyH'
 end select 

 call HiggsSignals_info 
 if(inputmethod=='subrout') then 
  if(allocated(theo))then
   if(debug) write(*,*) "HiggsBounds/HiggsSignals internal structure already initialized!"
  else
   if(debug)write(*,*)'doing other preliminary tasks...'      ; call flush(6)
   call setup_input_for_hs

   allocate(inputsub( 2 )) !(1)np(Hneut)>0 (2)np(Hplus)>0
   inputsub(1)%desc='HiggsBounds_neutral_input_*'    ; inputsub(1)%req=req(   0,   1)
   inputsub(2)%desc='HiggsBounds_charged_input'      ; inputsub(2)%req=req(   1,   0)
 
   do i=1,ubound(inputsub,dim=1)
    inputsub(i)%stat=0
   enddo
  endif
 endif 
 
 if(debug)write(*,*)'reading in Standard Model tables...'   ; call flush(6) 
 if(.not.allocated(BRSM)) call setup_BRSM 
 call setup_uncertainties
             
 if(debug)write(*,*)'reading in experimental data...'       ; call flush(6)
 call setup_observables

 if(debug)write(*,*)'sorting out processes to be checked...'; call flush(6)
 nanalys = size(obs)

 if(debug)write(*,*)'preparing output arrays...'            ; call flush(6)
 call setup_output_for_hs

 if(debug)write(*,*)'HiggsSignals has been initialized...'  ; call flush(6)

 just_after_run=.False.
  
 contains 
   !         |   np
   !         |Hneu Hcha 
   !         | ==0  ==0 
 function req(Hneu,Hcha)
  integer, intent(in) ::Hneu,Hcha
  integer :: req
  
  req=1
  if(np(Hneut)==0)  req= Hneu  * req
  if(np(Hplus)==0)  req= Hcha  * req

 end function req 

end subroutine initialize_HiggsSignals
!------------------------------------------------------------
subroutine HiggsSignals_neutral_input_MassUncertainty(dMh)
! Sets the theoretical mass uncertainty of the Higgs bosons.
!------------------------------------------------------------
 use usefulbits, only: theo,np,Hneut
 
 implicit none
 double precision,intent(in) :: dMh(np(Hneut))

 if(.not.allocated(theo))then
  stop 'subroutine HiggsSignals_initialize must be called first'
 endif
 if(np(Hneut).eq.0)then
  write(*,*)'subroutine HiggsSignal_neutral_input_MassUncertainty should'
  write(*,*)'only be called if np(Hneut)>0'
  stop 'error in subroutine HiggsSignal_neutral_input_MassUncertainty'
 endif
 
 theo(1)%particle(Hneut)%dM = dMh
 
end subroutine HiggsSignals_neutral_input_MassUncertainty
!------------------------------------------------------------
subroutine setup_uncertainties
!------------------------------------------------------------
 use usefulbits, only : file_id_common3
 use store_pathname_hs, only : pathname_HS
 use usefulbits_hs, only : delta_rate
 use io, only : read_matrix_from_file
 
 logical :: BRmodel, BRSM, XSmodel, XSSM

 call read_matrix_from_file(9,pathname_HS//"BRcov.in",delta_rate%BRcov, BRmodel)
 call read_matrix_from_file(9,pathname_HS//"BRcovSM.in",delta_rate%BRcovSM, BRSM)
 call read_matrix_from_file(5,pathname_HS//"XScov.in",delta_rate%CScov, XSmodel)
 call read_matrix_from_file(5,pathname_HS//"XScovSM.in",delta_rate%CScovSM, XSSM)

 if(BRmodel.and.BRSM) then
  delta_rate%BRcov_ok=.True.
  write(*,*) "Covariance matrix for relative branching ratio uncertainties read in successfully."
 else
  write(*,*) "Covariance matrix for relative branching ratio uncertainties not provided. Using default values."
 endif
 if(XSmodel.and.XSSM) then
  delta_rate%CScov_ok=.True.
  write(*,*) "Covariance matrix for relative cross section uncertainties read in successfully."
 else
  write(*,*) "Covariance matrix for relative cross section uncertainties not provided. Using default values."
 endif
    
end subroutine setup_uncertainties
!------------------------------------------------------------
subroutine setup_rate_uncertainties( dCS, dBR )
!------------------------------------------------------------
! Sets (relative) systematic uncertainties of the model for:
!  dCS(1) - singleH				dBR(1) - gamma gamma
!  dCS(2) - VBF					dBR(2) - W W
!  dCS(3) - HW					dBR(3) - Z Z
!  dCS(4) - HZ					dBR(4) - tau tau
!  dCS(5) - ttH					dBR(5) - b bbar
!------------------------------------------------------------
 use usefulbits_hs, only : delta_rate
 implicit none
 
 double precision, intent(in) :: dCS(5)
 double precision, intent(in) :: dBR(5)
 integer :: i
   
 delta_rate%dCS = dCS

 do i=lbound(dBR,dim=1),ubound(dBR,dim=1) 
  call setup_dbr(i,dBR(i))
 enddo
 
end subroutine setup_rate_uncertainties
!------------------------------------------------------------
subroutine setup_dbr(BRid, value)
!------------------------------------------------------------
 use usefulbits_hs, only : delta_rate

 integer,intent(in) :: BRid
 double precision, intent(in) :: value

 if(BRid.gt.0.and.BRid.lt.10) then
  delta_rate%dBR(BRid) = value
 else
  write(*,*) "Warning in setup_dbr: Unknown decay mode."
 endif
   
end subroutine setup_dbr
!------------------------------------------------------------
subroutine setup_correlations(corr_mu, corr_mh)
!------------------------------------------------------------
! With this subroutine the user may switch off/on correlations
! (default=on) by setting corr = 0/1.
!------------------------------------------------------------
 use usefulbits_hs, only : correlations_mu, correlations_mh
 implicit none
 
 integer, intent(in) :: corr_mu, corr_mh
 if(corr_mu.eq.0) then
  correlations_mu = .False.
  write(*,*) 'Correlations in signal strength observables are switched off.'
 elseif(corr_mu.eq.1) then
  correlations_mu = .True.
 else
  stop 'Error: Correlations must be switched on/off by an integer value of 0 or 1.' 
 endif
 if(corr_mh.eq.0) then
  correlations_mh = .False.
  write(*,*) 'Correlations in Higgs mass observables are switched off.'
 elseif(corr_mh.eq.1) then
  correlations_mh = .True.
 else
  stop 'Error: Correlations must be switched on/off by an integer value of 0 or 1.' 
 endif
end subroutine setup_correlations
!------------------------------------------------------------
subroutine setup_symmetricerrors(symm)
! Sets the measured rate uncertainties to either a symmetrical average
! of the upper and lower cyan band widths (symm==1) or else uses the 
! original (asymmetrical) errors.
!------------------------------------------------------------
 use usefulbits_hs, only : symmetricerrors
 implicit none

 integer, intent(in) :: symm
 if(symm.eq.1) then
  write(*,*) "Using averaged (symmetrical) experimental rate uncertainties."
  symmetricerrors = .True.
 else
  write(*,*) "Using original (asymmetrical) experimental rate uncertainties."
  symmetricerrors = .False.
 endif
  
end subroutine setup_symmetricerrors
!------------------------------------------------------------
subroutine setup_absolute_errors(absol)
! Treats the measured rate uncertainties as either absolute
! uncertainties (1) or relative (0). By default, they are
! treated as relative uncertainties.
!------------------------------------------------------------
 use usefulbits_hs, only : absolute_errors
 implicit none

 integer, intent(in) :: absol
 if(absol.eq.1) then
  write(*,*) "Using absolute experimental rate uncertainties."
  absolute_errors = .True.
 else
  write(*,*) "Using relative experimental rate uncertainties."
  absolute_errors = .False.
 endif
  
end subroutine setup_absolute_errors
!------------------------------------------------------------
subroutine setup_correlated_rate_uncertainties(corr)
!------------------------------------------------------------
use usefulbits_hs, only : delta_rate
integer, intent(in) :: corr

if(corr.eq.0) then
 delta_rate%usecov = .False.
 write(*,*) "Deactivated correlated CS and BR uncertainties. Using approximated maximum error."
elseif(corr.eq.1) then
 delta_rate%usecov = .True.
 write(*,*) "Activated correlated CS and BR uncertainties. Using them if covariance matrices are present."
else
 write(*,*) "Warning in subroutine setup_correlated_rate_uncertainties: Argument ",corr," is not equal to 0 or 1."
endif
end subroutine setup_correlated_rate_uncertainties
!------------------------------------------------------------
subroutine setup_SMweights(useweight)
! If set to 1 (true), HiggsSignals assumes the same signal decomposition
! (weights) as in the SM for the given model. This will enter the determination
! of the theoretical rate uncertainty.
!------------------------------------------------------------
 use usefulbits_hs, only : useSMweights
 implicit none

 integer, intent(in) :: useweight
 if(useweight.eq.1) then
  write(*,*) "Using SM weights for theoretical rate uncertainties of the model."
  useSMweights = .True.
 else
  write(*,*) "Using true model weights for theoretical rate uncertainties of the model."
  useSMweights = .False.
 endif
  
end subroutine setup_SMweights
!------------------------------------------------------------
subroutine setup_anticorrelations_in_mu(acorr)
! Allows for anti-correlations in the signal strength covariance
! matrix if there is a relative sign difference in two mu measurements
! (acorr==1) or else uses only correlations irrespective of the relative
! (acorr==0).
!------------------------------------------------------------
 use usefulbits_hs, only : anticorrmu
 implicit none

 integer, intent(in) :: acorr
 if(acorr.eq.1) then
  write(*,*) "Allow anti-correlated signal strength measurements."
  anticorrmu = .True.
 else
  write(*,*) "Prohibit anti-correlated signal strength measurements."
  anticorrmu = .False.
 endif
  
end subroutine setup_anticorrelations_in_mu
!------------------------------------------------------------
subroutine setup_anticorrelations_in_mh(acorr)
! Allows for anti-correlations in the mass covariance
! matrix if there is a relative sign difference in two mu measurements
! (acorr==1) or else uses only correlations irrespective of the relative
! (acorr==0).
!------------------------------------------------------------
 use usefulbits_hs, only : anticorrmh
 implicit none

 integer, intent(in) :: acorr
 if(acorr.eq.1) then
  write(*,*) "Allow anti-correlated mass measurements."
  anticorrmh = .True.
 else
  write(*,*) "Prohibit anti-correlated mass measurements."
  anticorrmh = .False.
 endif
  
end subroutine setup_anticorrelations_in_mh
!------------------------------------------------------------
subroutine setup_assignmentrange(range)
!------------------------------------------------------------
! This sets up the mass range (in standard deviations) in which
! the Higgs is forced to be assigned to the peak observables.
!------------------------------------------------------------
 use usefulbits_hs, only : assignmentrange,assignmentrange_massobs, pdf
 implicit none
 
 double precision, intent(in) :: range
 
 if(range.le.0.0D0) then
  write(*,*) "Error: Bad assignment range ",range
  write(*,*) "Keeping the value ",assignmentrange  
 else
  assignmentrange = range
  assignmentrange_massobs = range
 endif 

 if(assignmentrange.ne.1.0D0.and.pdf.eq.1) then
  write(*,*) "Note: For a box pdf, only 1s mass range is used to force the Higgs-to-peak assignment."
 endif
 
end subroutine setup_assignmentrange
!------------------------------------------------------------
subroutine setup_assignmentrange_massobservables(range)
!------------------------------------------------------------
! This sets up the mass range (in standard deviations) in which
! the Higgs is forced to be assigned to the peak observables.
!------------------------------------------------------------
 use usefulbits_hs, only : assignmentrange_massobs, pdf
 implicit none
 
 double precision, intent(in) :: range
 
 if(range.le.0.0D0) then
  write(*,*) "Error: Bad assignment range ",range
  write(*,*) "Keeping the value ",assignmentrange_massobs
 else
  assignmentrange_massobs = range
 endif 

 if(assignmentrange_massobs.ne.1.0D0.and.pdf.eq.1) then
  write(*,*) "Note: For a box pdf, only 1s mass range is used to force the Higgs-to-peak assignment."
 endif
 
end subroutine setup_assignmentrange_massobservables
!------------------------------------------------------------
subroutine setup_nparam(Np)
!------------------------------------------------------------
 use usefulbits_hs, only : Nparam
 implicit none 
 integer, intent(in) :: Np
 Nparam = Np
end subroutine setup_nparam
!------------------------------------------------------------
subroutine setup_Higgs_to_peaks_assignment_iterations(iter)
! Sets the number of iterations for the Higgs-to-peak-assignment.
!------------------------------------------------------------
 use usefulbits_hs, only : iterations
 implicit none
 integer, intent(in) :: iter 
 iterations = iter 
 
end subroutine setup_Higgs_to_peaks_assignment_iterations
!------------------------------------------------------------
subroutine setup_mcmethod_dm_theory(mode)
 use mc_chisq, only : mc_mode
 implicit none
 integer, intent(in) :: mode
 character(LEN=14) :: mode_desc(2) = (/'mass variation','convolution   '/)

 if(mode.eq.1.or.mode.eq.2) then
  mc_mode = mode
  write(*,'(1X,A,A)') 'The mass-centered chi^2 method will treat the Higgs',&
& ' boson mass theory uncertainty by '//trim(mode_desc(mode))//'.'
 else
  stop 'Error in subroutine setup_mcmethod_dm_theory: Unknown mode (1 or 2 possible)!'
 endif 
end subroutine setup_mcmethod_dm_theory
!------------------------------------------------------------
subroutine setup_sm_test(int_SMtest,epsilon)
! With this subroutine the user may switch off the SM likeness test
! (default=on) or change the maximal deviation epsilon (default=5.0D-2)
!------------------------------------------------------------
 use usefulbits_hs, only : useSMtest, eps
 implicit none
 
 integer, intent(in) :: int_SMtest
 double precision, intent(in) :: epsilon
 
 if(int_SMtest.eq.0) then
  useSMtest = .False.
  write(*,*) 'SM likeness test has been switched off.'
 elseif(int_SMtest.eq.1) then
  useSMtest = .True.
  write(*,*) 'SM likeness test has been switched on.'  
 else
  stop 'Error: SM test must be switched on/off by an integer value of 0 or 1.' 
 endif
 eps = epsilon
end subroutine setup_sm_test
!------------------------------------------------------------
subroutine setup_thu_observables(thuobs)
 use usefulbits_hs, only : THU_included
 integer, intent(in) :: thuobs

 if(thuobs.eq.0) then
  THU_included = .False.
  write(*,*) 'Observables are assumed to NOT include theory errors.'
 else
  THU_included = .True.
  write(*,*) 'Observables are assumed to include theory errors.'  
 endif
  
end subroutine setup_thu_observables
!------------------------------------------------------------
subroutine setup_output_level(level)
! Controls the level of information output:
! 0 : silent mode
! 1 : screen output for each analysis with its peak/mass-centered observables and
!     their respective values predicted by the model
! 2 : screen output of detailed information on each analysis with its 
!     peak/mass-centered observables
! 3 : creates the files peak_information.txt and peak_massesandrates.txt
!------------------------------------------------------------
 use usefulbits_hs, only : output_level, additional_output
 implicit none
 integer, intent(in) :: level

 if(level.eq.0.or.level.eq.1.or.level.eq.2.or.level.eq.3) then
  output_level = level
 else
  stop 'Error in subroutine setup_output_level: level not equal to 0,1,2 or 3.' 
 endif
 if(level.eq.3) additional_output = .True.
 
end subroutine setup_output_level 
!------------------------------------------------------------
subroutine setup_pdf(pdf_in)
! Sets the probability density function for the Higgs mass uncertainty parametrization:
! 1 : box-shaped pdf
! 2 : Gaussian pdf
! 3 : box-shaped theory error + Gaussian experimental pdf
!------------------------------------------------------------
 use usefulbits_hs, only : pdf, assignmentrange

 implicit none
 integer, intent(in) :: pdf_in
 character(LEN=13) :: pdf_desc(3) = (/'box         ','Gaussian    ','box+Gaussian'/)
 
 pdf=pdf_in
 if((pdf.eq.1).or.(pdf.eq.2).or.(pdf.eq.3)) then
  write(*,'(1X,A,A,1I1,A)') 'Use a '//trim(pdf_desc(pdf))//' probability density function ',&
&  'for the Higgs mass(es) (pdf=',pdf,')'
 endif
   
 if(assignmentrange.ne.1.0D0.and.pdf.eq.1) then
  write(*,*) "Note: For a box pdf, only 1s mass range is used to force the Higgs-to-peak assignment."
 endif
   
end subroutine setup_pdf
!------------------------------------------------------------
!subroutine assign_toyvalues_to_observables(ii, peakindex, npeaks, mu_obs, mh_obs)
!! Assigns toy values to the peak's mass and mu value for analysis ii.
!! ii           :: analysis number (entry in mutables)
!! peakindex    :: index of the peak of analysis ii
!! npeaks 	   :: number of peaks found in analysis ii
!! mu_obs       :: toy value for mu to be given to the peak with peakindex
!! mh_obs       :: toy value for mh to be given to the peak with peakindex
!------------------------------------------------------------
! use usefulbits_hs, only: obs, usetoys
! 
! integer, intent(in) :: ii, peakindex, npeaks
! double precision, intent(in) :: mh_obs, mu_obs
!  
! if(peakindex.gt.npeaks) then
!  stop 'Error in subroutine assign_toyvalues_to_observables: Observable does not exist!'
! endif
!   
! obs(ii)%table%npeaks = npeaks
! if(.not.allocated(obs(ii)%table%Toys_muobs)) allocate(obs(ii)%table%Toys_muobs(npeaks))
! if(.not.allocated(obs(ii)%table%Toys_mhobs)) allocate(obs(ii)%table%Toys_mhobs(npeaks)) 
!
! obs(ii)%table%Toys_muobs(peakindex) = mu_obs
! obs(ii)%table%Toys_mhobs(peakindex) = mh_obs
! 
! usetoys = .True.
! 
!end subroutine assign_toyvalues_to_observables
!------------------------------------------------------------
subroutine assign_toyvalues_to_peak(ID, mu_obs, mh_obs)
! Assigns toy values to the peak's mass and mu value to a peak observable.
! ID           :: observable ID
! mu_obs       :: toy value for mu to be given to the peak
! mh_obs       :: toy value for mh to be given to the peak
!
! n.B.: Do we also want to set mu uncertainties here?
!------------------------------------------------------------
 use usefulbits_hs, only: obs, usetoys
 implicit none

 integer, intent(in) :: ID
 double precision, intent(in) :: mh_obs, mu_obs
 integer :: pos, ii
 
 pos = -1
 do ii=lbound(obs,dim=1),ubound(obs,dim=1)
  if(obs(ii)%id.eq.ID) then
   pos = ii
   exit   	
  endif
 enddo
  
 if(pos.ne.-1) then    
  obs(pos)%peak%mpeak = mh_obs
  obs(pos)%peak%mu = mu_obs
  usetoys = .True.
 else
  write(*,*) "WARNING in assign_toyvalues_to_peak: ID unknown."
 endif
  
end subroutine assign_toyvalues_to_peak
!------------------------------------------------------------
subroutine assign_modelefficiencies_to_peak(ID, Nc, eff_ratios)
! Assigns to each channel of the observable the efficiency in the model 
! w.r.t the SM efficiency (as a ratio!)
! 
! ID           :: observable ID
! Nc			:: number of channels
! eff_ratios    :: array of length (Number of channels) giving the efficiency ratios
!
! Note: You can first employ the subroutine get_peak_channels (io module) to obtain
!       the relevant channel information of the observable.
!------------------------------------------------------------
 use usefulbits_hs, only: obs
 implicit none

 integer, intent(in) :: ID, Nc
 double precision, dimension(Nc), intent(in) :: eff_ratios
 integer :: pos, ii
 
 pos = -1
 do ii=lbound(obs,dim=1),ubound(obs,dim=1)
  if(obs(ii)%id.eq.ID) then
   pos = ii
   exit   	
  endif
 enddo
  
 if(pos.ne.-1) then    
  if(size(eff_ratios,dim=1).ne.obs(pos)%table%Nc) then
   write(*,*) "WARNING in assign modelefficiencies_to_peak: Number of channels (",&
&  size(eff_ratios,dim=1),"!=",obs(pos)%table%Nc,"does not match for observable ID = ",ID
  else
   obs(pos)%table%channel_eff_ratios = eff_ratios
  endif 
 else
  write(*,*) "WARNING in assign_modelefficiencies_to_peak: ID unknown."
 endif
  
end subroutine assign_modelefficiencies_to_peak
!------------------------------------------------------------
subroutine assign_rate_uncertainty_scalefactor_to_peak(ID, scale_mu)
! Assigns a rate uncertainty scalefactor to the peak specified by ID.
! This scalefactor will only scale the experimental rate uncertainties.
! The theory rate uncertainties must be given manually via setup_rate_uncertainties.
!
! ID       :: observable ID of the peak observable
! scale_mu :: scale_mu by which the mu uncertainty is scaled
!------------------------------------------------------------
 use usefulbits_hs, only: obs, usescalefactor
 implicit none
 
 integer, intent(in) :: ID
 double precision, intent(in) :: scale_mu
 integer :: pos, ii
 
 pos = -1
 do ii=lbound(obs,dim=1),ubound(obs,dim=1)
  if(obs(ii)%id.eq.ID) then
   pos = ii
   exit   	
  endif
 enddo
  
 if(pos.ne.-1) then
  obs(pos)%peak%scale_mu = scale_mu
 else
  write(*,*) "WARNING in assign_uncertainty_scalefactors_to_peak: ID unknown."
 endif
 usescalefactor = .True.
 
end subroutine assign_rate_uncertainty_scalefactor_to_peak
!------------------------------------------------------------
subroutine run_HiggsSignals(mode, Chisq_mu, Chisq_mh, Chisq, nobs, Pvalue)
!------------------------------------------------------------
! This subroutine can be called by the user after HiggsSignals_initialize has been called.
! The input routines, where required, should be called before calling run_HiggsSignals.
! It takes theoretical predictions for a particular parameter point 
! in the model and calls subroutines which compare these predictions 
! to the experimental results.
! Arguments (output):
!   * mode = 1,2 or 3 for peak-centered, mass-centered chi^2 method or both, respectively. 
!   * Chisq_mu = total chi^2 contribution from signal strength measurements
!   * Chisq_mh = total chi^2 contribution from Higgs mass measurements
!   * Chisq = total chi^2 value for the combination of the considered Higgs signals
!   * nobs = total number of observables
!   * Pvalue = total chi^2 probability for the agreement between model and data,
!              assuming number of observables == number of degrees of freedom
!    (see manual for more precise definitions))
!------------------------------------------------------------
 use usefulbits, only : theo,inputsub,just_after_run, inputmethod, ndat
 use usefulbits_HS, only : HSres, runmode, output_level, usescalefactor, Nparam
 use channels, only : check_channels
 use theo_manip, only : complete_theo

#if defined(NAGf90Fortran)
 use F90_UNIX_IO, only : flush
#endif

 implicit none               
 integer,intent(in) :: mode 
 !----------------------------------------output
 integer,intent(out) ::           nobs
 double precision,intent(out) ::  Pvalue, Chisq, Chisq_mu, Chisq_mh
 !-------------------------------------internal
 integer :: n,i
 logical :: debug=.False.
 !---------------------------------------------
 
if(mode.eq.1) then
 runmode="peak"
else if(mode.eq.2) then
 runmode="mass"
else if(mode.eq.3) then
 runmode="both"
else 
 stop'Error in subroutine run_HiggsSignals: mode unknown'
endif

 if(.not.allocated(theo))then
  stop 'subroutine HiggsSignals_initialize must be called first'
 endif
 if(inputmethod.eq.'subrout') then
  do i=1,ubound(inputsub,dim=1)
   if(  inputsub(i)%req .ne. inputsub(i)%stat  )then
!    write(*,*) inputsub(i)%req, inputsub(i)%stat
!    write(*,*)'subroutine '//trim(adjustl(inputsub(i)%desc))
!    write(*,*)'should be called once and only once before each call to'
!    write(*,*)'subroutine run_HiggsSignals.'
!    stop'error in subroutine run_HiggsSignals'
   endif
! TS: Have to work on this bit to make it run simultaneously with HiggsBounds. Now,
!     commented out the =0 statement. HS thus has to be run before HB.   
  inputsub(i)%stat=0!now we have used this input, set back to zero   
  enddo
 endif 

 if(debug)write(*,*)'manipulating input...'                 ; call flush(6)

 call complete_theo       

 if(debug)write(*,*)'compare each model to the experimental data...' ; call flush(6)                  

 do n=1,ndat

  call evaluate_model(theo(n),HSres(n))       
      
  Pvalue  = HSres(n)%Pvalue
  Chisq   = HSres(n)%Chisq       
  Chisq_mu   = HSres(n)%Chisq_mu       
  Chisq_mh   = HSres(n)%Chisq_mh         
  nobs = HSres(n)%nobs

 if(output_level.ne.0) then    
  write(*,*)
  write(*,*) '#*************************************************************************#' 
  write(*,*) '#                         HIGGSSIGNALS RESULTS                            #'
  write(*,*) '#*************************************************************************#' 
  write(*,'(A55,F21.8)') 'chi^2 from signal strength peak observables = ',&
 &  HSres(n)%Chisq_peak_mu
  write(*,'(A55,F21.8)') 'chi^2 from Higgs mass peak observables = ',HSres(n)%Chisq_mh
  write(*,'(A55,F21.8)') 'chi^2 from mass-centered observables = ',HSres(n)%Chisq_mpred
  write(*,'(A55,F21.8)') 'chi^2 from signal strength (total) = ',HSres(n)%Chisq_mu
  write(*,'(A55,F21.8)') 'chi^2 (total) = ',HSres(n)%Chisq
  write(*,'(A55,I21)') 'Number of signal strength peak observables = ',&
 &  HSres(n)%nobs_peak_mu
  write(*,'(A55,I21)') 'Number of Higgs mass peak observables = ',HSres(n)%nobs_peak_mh
  write(*,'(A55,I21)') 'Number of mass-centered observables = ',HSres(n)%nobs_mpred
  write(*,'(A55,I21)') 'Number of observables (total) = ',HSres(n)%nobs
  write(*,'(A48,I3,A4,F21.8)') 'Probability (ndf =',HSres(n)%nobs-Nparam,') = ',HSres(n)%Pvalue
  write(*,*) '#*************************************************************************#' 
  write(*,*)
 endif

 enddo

 just_after_run=.True.
 usescalefactor=.False.
 

end subroutine run_HiggsSignals
!------------------------------------------------------------
subroutine evaluate_model( t , r )
!------------------------------------------------------------
! This subroutine evaluates the signal strength modifier for every Higgs boson and
! considered analysis. It fills a matrix neutHiggs(:,:) of type neutHiggs with dimensions
! (number(considered analyses),nH).
!------------------------------------------------------------
 use usefulbits, only : np,Hneut,Hplus,dataset,results, vsmall
 use usefulbits_hs, only : neutHiggses, nanalys, runmode, HSresults, cov, obs, analyses,&
 &						   cov_mhneut, iterations, deallocate_covariance_matrices, &
 &						   output_level, Nparam, nanalys
 use datatables, only : setup_tablelist, check_available_Higgses
 use pc_chisq
 use mc_chisq
 use all_chisq
 use numerics
 implicit none
 !--------------------------------------input      
 type(dataset), intent(in) :: t      
 !-------------------------------------output
 type(HSresults), intent(out) :: r

 integer :: ii, jj, iii, jjj
 
 double precision :: totchisq, muchisq, mhchisq, mpchisq, mpredchisq
 integer :: nobs, Nmu, Nmh, Nmpred
 character(LEN=100), allocatable :: assignmentgroups(:)
 integer, allocatable :: assignmentgroups_domH(:)
 integer, allocatable :: assignmentgroups_Higgs_comb(:,:)

 allocate(assignmentgroups(nanalys),assignmentgroups_domH(nanalys))
 allocate(assignmentgroups_Higgs_comb(nanalys,np(Hneut))) 

!---Initialize assignmentgroups arrays with default values
do ii=lbound(assignmentgroups_domH,dim=1),ubound(assignmentgroups_domH,dim=1)
 assignmentgroups_domH(ii) = 0
 assignmentgroups_Higgs_comb(ii,:) = 0
enddo    
  
!---First, evaluate the model predictions  
 allocate(neutHiggses(nanalys,np(Hneut)))
!-Loop over considered analyses  
 do ii=lbound(neutHiggses,dim=1),ubound(neutHiggses,dim=1)
!-Loop over the neutral Higgs bosons of the model 
  do jj=lbound(neutHiggses,dim=2),ubound(neutHiggses,dim=2)
!!   write(*,*) "hello evaluate model:", ii, jj
   call calc_mupred(jj, t, obs(ii)%table, neutHiggses(ii,jj))
  enddo
  if(.not.allocated(obs(ii)%Higgses)) allocate(obs(ii)%Higgses(np(Hneut)))
  obs(ii)%Higgses(:) = neutHiggses(ii,:)  
 enddo

!-Pass the observables and their predicted Higgs properties (obs%Higgses)
!-to the tablelist "analyses"
 call setup_tablelist

 select case(runmode)
 
 case('peak')
!-Peak-centered chisq method 
  jjj=0
  do ii=lbound(analyses,dim=1),ubound(analyses,dim=1)
   call deallocate_covariance_matrices
   call assign_Higgs_to_peaks(analyses(ii)%table, analyses(ii)%peaks,0)
   do iii=lbound(analyses(ii)%peaks,dim=1),ubound(analyses(ii)%peaks,dim=1)  
    if(analyses(ii)%table%mhchisq.eq.1.and.&
&      len(trim(adjustl(analyses(ii)%peaks(iii)%assignmentgroup))).ne.0) then
     jjj=jjj+1
     assignmentgroups(jjj)=analyses(ii)%peaks(iii)%assignmentgroup
     assignmentgroups_Higgs_comb(jjj,:)=analyses(ii)%peaks(iii)%Higgs_comb
     assignmentgroups_domH(jjj)=analyses(ii)%peaks(iii)%domH     
!!     write(*,*) "Found leader of group ",assignmentgroups(jjj)
!!     write(*,*) "ID ",analyses(ii)%peaks(iii)%id
!!     write(*,*) "with Higgs combination ",assignmentgroups_Higgs_comb(jjj,:)  
!!     write(*,*) "and dominant Higgs boson ",assignmentgroups_domH(jjj)          
    endif   
   enddo 
  enddo
  do ii=lbound(analyses,dim=1),ubound(analyses,dim=1)
   do iii=lbound(analyses(ii)%peaks,dim=1),ubound(analyses(ii)%peaks,dim=1)  
    if(analyses(ii)%table%mhchisq.eq.0.and.&
&     len(trim(adjustl(analyses(ii)%peaks(iii)%assignmentgroup))).ne.0) then
      !SELECT ASSIGNMENT GROUP FOLLOWERS
      do jjj=lbound(assignmentgroups,dim=1),ubound(assignmentgroups,dim=1)
       if(analyses(ii)%peaks(iii)%assignmentgroup.eq.assignmentgroups(jjj)) then
        !TAKE OVER THE HIGGS ASSIGNMENT OF THE LEADING PEAK
        analyses(ii)%peaks(iii)%Higgs_comb=assignmentgroups_Higgs_comb(jjj,:)
        analyses(ii)%peaks(iii)%domH=assignmentgroups_domH(jjj)
        if(assignmentgroups_domH(jjj).ne.0) then
         analyses(ii)%peaks(iii)%Higgs_assignment_forced=1
        endif 
        call evaluate_peak(analyses(ii)%peaks(iii),analyses(ii)%table)
       endif
      enddo
    endif
   enddo
  enddo     

! Do the iterative Higgs-to-peak-assignment here:
  call assign_Higgs_to_peaks_with_correlations(iterations)
  call calculate_total_pc_chisq(totchisq, muchisq, mhchisq, nobs, Nmu, Nmh)

  if(output_level.eq.1) call print_peakinformation
  if(output_level.eq.2) call print_peakinformation_essentials
  if(output_level.eq.3) then
   call print_peaks_to_file
   call print_peaks_signal_rates_to_file
  endif 

  call add_peaks_to_HSresults(r)
    
  r%Chisq=totchisq
  r%Chisq_peak_mu = muchisq
  r%Chisq_mpred = 0.0D0
  r%Chisq_mu=muchisq
  r%Chisq_mh=mhchisq  
  r%nobs_mpred=0
  r%nobs_peak_mu=Nmu
  r%nobs_peak_mh=Nmh
  r%nanalysis=size(analyses)
  r%nobs=nobs
  if(r%Chisq.gt.vsmall.and.(r%nobs-Nparam).gt.0) then
   r%Pvalue=1 - gammp(dble(r%nobs-Nparam)/2,r%Chisq/2)
  endif
  
 case('mass')  
  do ii=lbound(analyses,dim=1),ubound(analyses,dim=1)
   call fill_mp_obs(ii)
  enddo
  if(mc_mode.eq.1) call mass_variation_by_theory_uncertainty
  call create_covariance_matrix_mp
  call calculate_mpred_chisq(mpchisq, nobs)

  if(output_level.eq.1) call print_mc_observables
  if(output_level.eq.2) call print_mc_observables_essentials
  if(output_level.eq.3) then
   call print_mc_tables_to_file
   call print_mc_observables_to_file   
  endif
  
  r%Chisq=mpchisq
  r%Chisq_peak_mu = 0.0D0
  r%Chisq_mpred = mpchisq 
  r%Chisq_mu=mpchisq
  r%Chisq_mh=0.0D0
  r%nobs_mpred=nobs
  r%nobs_peak_mu=0
  r%nobs_peak_mh=0
  r%nanalysis=size(analyses)  
  r%nobs=nobs    
  if(r%Chisq.gt.vsmall.and.(r%nobs-Nparam).gt.0) then
   r%Pvalue=1 - gammp(dble(r%nobs-Nparam)/2,r%Chisq/2)
  endif

 case('both')
 jjj=0
  do ii=lbound(analyses,dim=1),ubound(analyses,dim=1)
   call deallocate_covariance_matrices
   call assign_Higgs_to_peaks(analyses(ii)%table, analyses(ii)%peaks,0)
   do iii=lbound(analyses(ii)%peaks,dim=1),ubound(analyses(ii)%peaks,dim=1)  
    if(analyses(ii)%table%mhchisq.eq.1.and.&
&      len(trim(analyses(ii)%peaks(iii)%assignmentgroup)).ne.0) then
     jjj=jjj+1
     assignmentgroups(jjj)=analyses(ii)%peaks(iii)%assignmentgroup
     assignmentgroups_Higgs_comb(jjj,:)=analyses(ii)%peaks(iii)%Higgs_comb
     assignmentgroups_domH(jjj)=analyses(ii)%peaks(iii)%domH     
    endif   
   enddo    
  enddo
  do ii=lbound(analyses,dim=1),ubound(analyses,dim=1)
   do iii=lbound(analyses(ii)%peaks,dim=1),ubound(analyses(ii)%peaks,dim=1)  
    if(analyses(ii)%table%mhchisq.eq.0.and.&
&     len(trim(analyses(ii)%peaks(iii)%assignmentgroup)).ne.0) then
      do jjj=lbound(assignmentgroups,dim=1),ubound(assignmentgroups,dim=1)
       if(analyses(ii)%peaks(iii)%assignmentgroup.eq.assignmentgroups(jjj)) then
        !TAKE OVER THE HIGGS ASSIGNMENT OF THE LEADING PEAK
        analyses(ii)%peaks(iii)%Higgs_comb=assignmentgroups_Higgs_comb(jjj,:)
        analyses(ii)%peaks(iii)%domH=assignmentgroups_domH(jjj)
        if(assignmentgroups_domH(jjj).ne.0) then
         analyses(ii)%peaks(iii)%Higgs_assignment_forced=1
        endif 
        ! TODO: Need to evaluate everything else here!
        call evaluate_peak(analyses(ii)%peaks(iii),analyses(ii)%table)
       endif
      enddo
    endif
   enddo
  enddo     
  
  call assign_Higgs_to_peaks_with_correlations(iterations) 
  
  do ii=lbound(analyses,dim=1),ubound(analyses,dim=1)
   call check_available_Higgses(ii)
   call fill_mp_obs(ii)
  enddo  
  if(mc_mode.eq.1) call mass_variation_by_theory_uncertainty
 
  call calculate_total_chisq(totchisq, muchisq, mhchisq, mpredchisq, nobs, Nmu, Nmh, Nmpred)
 
 !Have to write a new print method
  if(output_level.eq.1)  call print_all_observables
  if(output_level.eq.2) call print_peakinformation_essentials
  if(output_level.eq.3) then
   call print_peaks_to_file
   call print_peaks_signal_rates_to_file
  endif 

  call add_peaks_to_HSresults(r)
 
  r%Chisq=totchisq
  r%Chisq_peak_mu = muchisq
  r%Chisq_mpred = mpredchisq 
  r%Chisq_mu=muchisq + mpredchisq
  r%Chisq_mh=mhchisq  
  r%nobs_mpred=Nmpred
  r%nobs_peak_mu=Nmu
  r%nobs_peak_mh=Nmh
  r%nanalysis=size(analyses)
  r%nobs=nobs
  if(r%Chisq.gt.vsmall.and.(r%nobs-Nparam).gt.0) then
   r%Pvalue=1 - gammp(dble(r%nobs-Nparam)/2,r%Chisq/2)
  endif
  
 case default
  stop "Error in subroutine evaluate_model: Please specify runmode!"
  
 end select

 deallocate(neutHiggses)
 deallocate(assignmentgroups, assignmentgroups_domH, assignmentgroups_Higgs_comb)   
end subroutine evaluate_model
!------------------------------------------------------------
subroutine calc_mupred( j, t, mutab, Higgs )
! Calculates the model-predicted signal strength modifier
!------------------------------------------------------------
 use usefulbits, only : dataset, div, vsmall
 use usefulbits_HS, only : neutHiggs, mutable, useSMtest, eps
 implicit none
 
 integer, intent(in) :: j					! Higgs index
 type(dataset), intent(in) :: t
 type(mutable), intent(inout) :: mutab
 type(neutHiggs), intent(inout) :: Higgs

 integer :: i
 double precision :: c, dcbyc
 integer :: testSMratios 
 logical :: correct_properties

 Higgs%m = t%particle(mutab%particle_x)%M(j)
 Higgs%dm = t%particle(mutab%particle_x)%dM(j)
 Higgs%id = j
  
 call get_channelrates( j, t, mutab )

 correct_properties=.True.

!--Evaluate the predicted signal strength modifier c of the model
 c=0. 
 do i=1,mutab%Nc
!----use a weighted average of the channel rate ratios     
  c=c+mutab%channel_w(i,j)*mutab%channel_mu(i,j)
 enddo

!--Evaluate the deviation of each channel rate ratio to the signal
!--strength modifier c and test SM likeness criterium, if this is
!--activated.
 testSMratios= 1  !passes the SM-like ratios test 
 do i=1,mutab%Nc
  dcbyc=div((mutab%channel_mu(i,j)-c),c,0.0D0,1.0D9)
  if(dcbyc*mutab%channel_w(i,j).gt.eps.and.useSMtest) then
   testSMratios= -1  !fails the SM-like ratios test
  endif     
 enddo

 if(testSMratios.lt.0) correct_properties=.False.
  
 if(correct_properties) then
  Higgs%mu=c
 else
  Higgs%mu=0.0D0
 endif
  
end subroutine calc_mupred
!------------------------------------------------------------
subroutine get_channelrates( j, t, mutab )
! This subroutine assignes the rates, weights and systematic rate uncertainty of
! the Higgs boson (j) for the channels considered by the analysis (mutab).
!
! WARNING: if normalize_rates_to_reference_position is true -> Still in TESTING PHASE!
! The rates are normalized w.r.t. a reference rate at the (peak) mass position.
! This does not work with the mass-centered chi^2 method.
! Also, theoretical mass uncertainties are problematic!
!------------------------------------------------------------
 use usefulbits, only : dataset, div, small
 use usefulbits_HS, only : neutHiggs, mutable, delta_rate, normalize_rates_to_reference_position
 use theory_XS_SM_functions
 use theory_BRfunctions
 
 integer, intent(in) :: j
 type(dataset), intent(in) :: t
 type(mutable), intent(inout) :: mutab


 integer :: i, id, p, d
 integer :: ii, id1, id2, p1, p2, d1, d2
 double precision :: rate, SMrate, modelrate, drsq_SM, drsq, dBR, dBRSM
 
!!NEW:
 double precision :: rate_SMref,refmass

 if(size(mutab%mass,dim=1).eq.1) then
  refmass = mutab%mass(1)
 else
!  write(*,*) "mutab%id", mutab%id, "Mass measurements: ",size(mutab%mass,dim=1)
!  write(*,*) "mutab%particle_x = ", mutab%particle_x, " j= ", j
  refmass = t%particle(mutab%particle_x)%M(j)
 endif

!! write(*,*) "hello: ", mutab%id

 do i=1,mutab%Nc
  id = mutab%channel_id(i)
  p = int((id-modulo(id,10))/dble(10))
  d = modulo(id,10)

!--Do the production rate for the relevant experiment and cms-energy 
  if(mutab%collider.eq.'LHC') then
   if(abs(mutab%energy-7.0D0).le.small) then
    if(p.eq.1) then 
     rate=t%lhc7%XS_hj_ratio(j)
     SMrate=t%lhc7%XS_H_SM(j)
     rate_SMref=XS_lhc7_gg_H_SM(refmass)
     mutab%channel_description(i,1)='singleH'
    else if(p.eq.2) then
     rate=t%lhc7%XS_vbf_ratio(j)
     SMrate=t%lhc7%XS_vbf_SM(j)
     rate_SMref=XS_lhc7_vbf_SM(refmass)
     mutab%channel_description(i,1)='VBF'     
    else if(p.eq.3) then
     rate=t%lhc7%XS_hjW_ratio(j)
     SMrate=t%lhc7%XS_HW_SM(j) 
     rate_SMref=XS_lhc7_HW_SM(refmass)
     mutab%channel_description(i,1)='HW'     
    else if(p.eq.4) then
     rate=t%lhc7%XS_hjZ_ratio(j)  
     SMrate=t%lhc7%XS_HZ_SM(j)
     rate_SMref=XS_lhc7_HZ_SM(refmass)
     mutab%channel_description(i,1)='HZ'       
    else if(p.eq.5) then
     rate=t%lhc7%XS_tthj_ratio(j)
     SMrate=t%lhc7%XS_ttH_SM(j)
     rate_SMref=XS_lhc7_ttH_SM(refmass)
     mutab%channel_description(i,1)='ttH'     
    else if(p.eq.0) then
     rate=1.0D0
     SMrate=1.0D0
     rate_SMref=1.0D0
     mutab%channel_description(i,1)='none'
    endif 
   else if(abs(mutab%energy-8.0D0).le.small) then
    if(p.eq.1) then 
     rate=t%lhc8%XS_hj_ratio(j)
     SMrate=t%lhc8%XS_H_SM(j)
     rate_SMref=XS_lhc8_gg_H_SM(refmass)     
     mutab%channel_description(i,1)='singleH'
    else if(p.eq.2) then
     rate=t%lhc8%XS_vbf_ratio(j)
     SMrate=t%lhc8%XS_vbf_SM(j)
     rate_SMref=XS_lhc8_vbf_SM(refmass)     
     mutab%channel_description(i,1)='VBF'     
    else if(p.eq.3) then
     rate=t%lhc8%XS_hjW_ratio(j)
     SMrate=t%lhc8%XS_HW_SM(j) 
     rate_SMref=XS_lhc8_HW_SM(refmass)     
     mutab%channel_description(i,1)='HW'     
    else if(p.eq.4) then
     rate=t%lhc8%XS_hjZ_ratio(j)  
     SMrate=t%lhc8%XS_HZ_SM(j)
     rate_SMref=XS_lhc8_HZ_SM(refmass)     
     mutab%channel_description(i,1)='HZ'       
    else if(p.eq.5) then
     rate=t%lhc8%XS_tthj_ratio(j)
     SMrate=t%lhc8%XS_ttH_SM(j)
     rate_SMref=XS_lhc8_ttH_SM(refmass)     
     mutab%channel_description(i,1)='ttH' 
    else if(p.eq.0) then
     rate=1.0D0
     SMrate=1.0D0
     rate_SMref=1.0D0     
     mutab%channel_description(i,1)='none'         
    endif  
   endif 
  else if(mutab%collider.eq.'TEV') then
    if(p.eq.1) then 
     rate=t%tev%XS_hj_ratio(j)
     SMrate=t%tev%XS_H_SM(j)
     rate_SMref=XS_tev_gg_H_SM(refmass)     
     mutab%channel_description(i,1)='singleH'
    else if(p.eq.2) then
     rate=t%tev%XS_vbf_ratio(j)
     SMrate=t%tev%XS_vbf_SM(j)
     rate_SMref=XS_tev_vbf_SM(refmass)     
     mutab%channel_description(i,1)='VBF'     
    else if(p.eq.3) then
     rate=t%tev%XS_hjW_ratio(j)
     SMrate=t%tev%XS_HW_SM(j) 
     rate_SMref=XS_tev_HW_SM(refmass)
     mutab%channel_description(i,1)='HW'     
    else if(p.eq.4) then
     rate=t%tev%XS_hjZ_ratio(j)  
     SMrate=t%tev%XS_HZ_SM(j)
     rate_SMref=XS_tev_HZ_SM(refmass)
     mutab%channel_description(i,1)='HZ'       
    else if(p.eq.5) then
     rate=t%tev%XS_tthj_ratio(j)
     SMrate=t%tev%XS_ttH_SM(j)
     rate_SMref=XS_tev_ttH_SM(refmass)
     mutab%channel_description(i,1)='ttH' 
    else if(p.eq.0) then
     rate=1.0D0
     SMrate=1.0D0
     rate_SMref=1.0D0     
     mutab%channel_description(i,1)='none'         
    endif       
  else if(mutab%collider.eq.'ILC') then
!--n.B.: As a first attempt, we use the LHC8 normalized cross sections for ZH, VBF, ttH.
!        In order to do this properly, a separate input for the ILC cross sections
!        has to be provided! It works only for single production mode observables (no
!        correct weighting of channels included!)Then, at least in the effective coupling
!        approximation, there is no difference to a full implementation.
!        The theoretical uncertainty of the ILC production modes will are defined in
!        usefulbits_HS.f90.
    if(p.eq.1.or.p.eq.2) then 
     write(*,*) 'Warning: Unknown ILC production mode (',p,') in table ',mutab%id
     rate=0.0D0
     SMrate=1.0D0
     rate_SMref=1.0D0     
     mutab%channel_description(i,1)='unknown'
    else if(p.eq.3) then
     rate=t%lhc8%XS_hjW_ratio(j)
     SMrate=t%lhc8%XS_HW_SM(j)
     rate_SMref=XS_lhc8_HW_SM(refmass)          
     mutab%channel_description(i,1)='WBF'     
    else if(p.eq.4) then
     rate=t%lhc8%XS_hjZ_ratio(j)  
     SMrate=t%lhc8%XS_HZ_SM(j)
     rate_SMref=XS_lhc8_HZ_SM(refmass)     
     mutab%channel_description(i,1)='HZ'       
    else if(p.eq.5) then
     rate=t%lhc8%XS_tthj_ratio(j)
     SMrate=t%lhc8%XS_ttH_SM(j)
     rate_SMref=XS_lhc8_ttH_SM(refmass)
     mutab%channel_description(i,1)='ttH' 
    else if(p.eq.0) then
     rate=1.0D0
     SMrate=1.0D0
     rate_SMref=1.0D0     
     mutab%channel_description(i,1)='none'         
    endif           
   endif
!--Multiply now by the decay rate
  if(d.eq.1) then
   rate=rate*div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0)
   SMrate=SMrate*t%BR_Hgaga_SM(j)
   rate_SMref = rate_SMref*BRSM_Hgaga(refmass)
   mutab%channel_description(i,2)='gammagamma'   
  else if(d.eq.2) then
   rate=rate*div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0)   
   SMrate=SMrate*t%BR_HWW_SM(j)
   rate_SMref = rate_SMref*BRSM_HWW(refmass)
   mutab%channel_description(i,2)='WW'   
  else if(d.eq.3) then
   rate=rate*div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0)
   SMrate=SMrate*t%BR_HZZ_SM(j)
   rate_SMref = rate_SMref*BRSM_HZZ(refmass)   
   mutab%channel_description(i,2)='ZZ'
  else if(d.eq.4) then
   rate=rate*div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0)
   SMrate=SMrate*t%BR_Htautau_SM(j)
   rate_SMref = rate_SMref*BRSM_Htautau(refmass)
   mutab%channel_description(i,2)='tautau'
  else if(d.eq.5) then
   rate=rate*div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0)
   SMrate=SMrate*t%BR_Hbb_SM(j)
   rate_SMref = rate_SMref*BRSM_Hbb(refmass)
   mutab%channel_description(i,2)='bb'  
  else if(d.eq.6) then
   rate=rate*div(t%BR_hjZga(j),t%BR_HZga_SM(j),0.0D0,1.0D0)
   SMrate=SMrate*t%BR_HZga_SM(j)
   rate_SMref = rate_SMref*BRSM_HZga(refmass)   
   mutab%channel_description(i,2)='Zgamma'
  else if(d.eq.7) then
   rate=rate*div(t%BR_hjcc(j),t%BR_Hcc_SM(j),0.0D0,1.0D0)
   SMrate=SMrate*t%BR_Hcc_SM(j)
   rate_SMref = rate_SMref*BRSM_Hcc(refmass)
   mutab%channel_description(i,2)='cc'
  else if(d.eq.8) then
   rate=rate*div(t%BR_hjmumu(j),t%BR_Hmumu_SM(j),0.0D0,1.0D0)
   SMrate=SMrate*t%BR_Hmumu_SM(j)
   rate_SMref = rate_SMref*BRSM_Hmumu(refmass)
   mutab%channel_description(i,2)='mumu'
  else if(d.eq.9) then
   rate=rate*div(t%BR_hjgg(j),t%BR_Hgg_SM(j),0.0D0,1.0D0)
   SMrate=SMrate*t%BR_Hgg_SM(j)
   rate_SMref = rate_SMref*BRSM_Hgg(refmass)
   mutab%channel_description(i,2)='gg'
  else if(d.eq.0) then
   rate=rate*1.0D0
   SMrate=SMrate*1.0D0
   rate_SMref = rate_SMref*1.0D0
   mutab%channel_description(i,2)='none'
  endif
  
 if(normalize_rates_to_reference_position) then 
!! THIS IS STILL IN TESTING PHASE !!  
  mutab%channel_mu(i,j)=rate*SMrate/(rate_SMref)
 else
  mutab%channel_mu(i,j)=rate  !! OLD WAY 
 endif
  
  mutab%channel_w(i,j)=mutab%channel_eff(i)*SMrate 
!  mutab%channel_w_corrected_eff(i,j)=mutab%channel_eff_ratios(i)*mutab%channel_eff(i)*SMrate   
 enddo

 SMrate=sum(mutab%channel_w(:,j))
! modelrate=sum(mutab%channel_w_corrected_eff(:,j))
 
 do i=1,mutab%Nc
  mutab%channel_w(i,j)=div(mutab%channel_w(i,j),SMrate,0.0D0,1.0D9)
!  mutab%channel_w_corrected_eff(i,j)=div(mutab%channel_w_corrected_eff(i,j),modelrate,0.0D0,1.0D9)
 enddo
 
! (TS 30/10/2013):
! write(*,*) "get_channelrates (mu, w, weff):" 
! write(*,*) mutab%channel_mu
! write(*,*)  mutab%channel_w
! write(*,*) mutab%channel_eff_ratios
 do i=1,mutab%Nc
  mutab%channel_w_corrected_eff(i,j)=mutab%channel_eff_ratios(i)*mutab%channel_w(i,j)
! n.b.: model weights are not normalized to 1!
 enddo
 
! write(*,*) j,mutab%id, "SM         = ", mutab%channel_w(:,j)
! write(*,*) j,mutab%id, "SM effcorr = ",mutab%channel_w_corrected_eff(:,j) 
 
 do i=1,mutab%Nc
  drsq_SM = 0.0D0
  drsq = 0.0D0

  id1 = mutab%channel_id(i)
  p1 = int((id1-modulo(id1,10))/dble(10))
  d1 = modulo(id1,10)
  if(mutab%collider.ne.'ILC') then
   do ii=1,mutab%Nc 
    id2 = mutab%channel_id(ii)
    p2 = int((id2-modulo(id2,10))/dble(10))
    d2 = modulo(id2,10)
    if(p1.eq.p2.and.p1.ne.0) then
     if(delta_rate%CScov_ok.and.delta_rate%usecov) then
      drsq=drsq+delta_rate%CScov(p1,p1)*mutab%channel_w_corrected_eff(i,j)*mutab%channel_w_corrected_eff(ii,j)
      drsq_SM=drsq_SM+delta_rate%CScovSM(p1,p1)*mutab%channel_w(i,j)*mutab%channel_w(ii,j)     
     else
      drsq=drsq+delta_rate%dCS(p1)**2*mutab%channel_w_corrected_eff(i,j)*mutab%channel_w_corrected_eff(ii,j)
      drsq_SM=drsq_SM+delta_rate%dCS_SM(p1)**2*mutab%channel_w(i,j)*mutab%channel_w(ii,j)
     endif
    endif 
    if(d1.eq.d2.and.d1.ne.0) then
     if(delta_rate%BRcov_ok.and.delta_rate%usecov) then
      dBRSM = delta_rate%BRcovSM(d1,d1)
      dBR = delta_rate%BRcov(d1,d1)
     else
      dBRSM = delta_rate%dBR_SM(d1)**2
      dBR = delta_rate%dBR(d1)**2
     endif 
     drsq=drsq+dBR*mutab%channel_w_corrected_eff(i,j)*mutab%channel_w_corrected_eff(ii,j)
     drsq_SM=drsq_SM+dBRSM*mutab%channel_w(i,j)*mutab%channel_w(ii,j)   
    endif 
   enddo 
  endif 
  mutab%channel_syst(i,j)=sqrt(drsq)
  mutab%channel_systSM(i,j)=sqrt(drsq_SM)
 enddo
   
end subroutine get_channelrates
!------------------------------------------------------------
subroutine get_Rvalues(ii,collider,R_H_WW, R_H_ZZ, R_H_gaga, R_H_tautau, R_H_bb, R_VH_bb)
! Returns SM normalized signal rates of some relevant channels (w/o efficiencies) 
! for Higgs boson "ii" for a specific collider (see subroutine get_rates).
!------------------------------------------------------------
! use usefulbits, only : theo, np,Hneut
! use usefulbits_HS, only : mutable

 integer, intent(in) :: ii, collider
 double precision, intent(out) :: R_H_WW, R_H_ZZ, R_H_gaga, R_H_tautau, R_H_bb, R_VH_bb
! type(mutable) :: dummytable
! integer :: i
 
 call get_rates(ii,collider,5,(/ 12, 22, 32, 42, 52 /),R_H_WW)
 call get_rates(ii,collider,5,(/ 13, 23, 33, 43, 53 /),R_H_ZZ)
 call get_rates(ii,collider,5,(/ 11, 21, 31, 41, 51 /),R_H_gaga)
 call get_rates(ii,collider,5,(/ 14, 24, 34, 44, 54 /),R_H_tautau)
 call get_rates(ii,collider,5,(/ 15, 25, 35, 45, 55 /),R_H_bb)
 call get_rates(ii,collider,2,(/ 35, 45 /),R_VH_bb)

end subroutine get_Rvalues
!************************************************************
subroutine get_rates(ii,collider,Nchannels,IDchannels,rate)
! Returns SM normalized signal rates (w/o efficiencies) for Higgs boson "ii" and collider
! experiment "collider"(=1,2,3 for TEV, LHC7, LHC8). "Nchannels" gives the total number
! and IDchannels the two-digit ID of the subchannels, which should be included in the rates.
! IDchannels is an array of size(Nchannels).
!------------------------------------------------------------
 use usefulbits, only : theo, np,Hneut
 use usefulbits_HS, only : mutable

 integer, intent(in) :: ii, collider, Nchannels
 integer, dimension(Nchannels), intent(in) :: IDchannels
 double precision, intent(out) :: rate
!-Internal 
 type(mutable) :: dummytable
 integer :: i
 
!-Initialize a dummy mutable in order to run get_channelrates for the channels we want.  
 if(collider.eq.1) then
  dummytable%collider = 'TEV'
 else if(collider.eq.2) then
  dummytable%collider = 'LHC'
  dummytable%energy = 7.0D0
 else if(collider.eq.3) then
  dummytable%collider = 'LHC'
  dummytable%energy = 8.0D0
 else
  write(*,*) 'WARNING: collider experiment for get_rates unknown.'
  continue
 endif 

 dummytable%id = 999999
 dummytable%particle_x = 1
 dummytable%Nc=Nchannels
 allocate(dummytable%mass(10))
 allocate(dummytable%channel_id(Nchannels))
 allocate(dummytable%channel_eff(Nchannels))
 allocate(dummytable%channel_eff_ratios(Nchannels))
!-Set all efficiencies equal: 
 dummytable%channel_eff = 1.0D0
 dummytable%channel_eff_ratios = 1.0D0
 allocate(dummytable%channel_description(Nchannels,2))
 allocate(dummytable%channel_w(Nchannels,np(Hneut)))
 allocate(dummytable%channel_w_corrected_eff(Nchannels,np(Hneut)))
 allocate(dummytable%channel_systSM(Nchannels,np(Hneut))) 
 allocate(dummytable%channel_syst(Nchannels,np(Hneut))) 
 allocate(dummytable%channel_mu(Nchannels,np(Hneut)))
 dummytable%channel_id = IDchannels

 call get_channelrates(ii, theo(1), dummytable)
 rate=0.0D0
 do i=lbound(dummytable%channel_mu,dim=1),ubound(dummytable%channel_mu,dim=1)
  rate = rate + dummytable%channel_mu(i,ii)*dummytable%channel_w(i,ii)
 enddo

 deallocate(dummytable%channel_id,dummytable%channel_eff,dummytable%channel_description,&
&           dummytable%channel_w,dummytable%channel_systSM,dummytable%channel_syst,		&
&           dummytable%channel_mu,dummytable%channel_eff_ratios, &
&           dummytable%channel_w_corrected_eff,dummytable%mass)

end subroutine get_rates
!------------------------------------------------------------
subroutine get_Pvalue(nparam, Pvalue)
! Calculates the Chi^2 probability for the total Chi^2 value
! and the number of degrees of freedom given by the
! number of observables - nparam
!------------------------------------------------------------
 use usefulbits, only : vsmall
 use usefulbits_hs, only: HSres
 use numerics 
 implicit none
 integer, intent(in) :: nparam
 double precision, intent(out) :: Pvalue
 
 if(allocated(HSres)) then
  if(HSres(1)%Chisq.gt.vsmall.and.(HSres(1)%nobs-nparam).gt.0) then
   HSres(1)%Pvalue = 1 - gammp(dble(HSres(1)%nobs-nparam)/2,HSres(1)%Chisq/2)
  endif
 else
  write(*,*) "Warning: subroutine get_Pvalue should be called after run_HiggsSignals." 
 endif
 
 Pvalue = HSres(1)%Pvalue

end subroutine get_Pvalue
!------------------------------------------------------------
subroutine get_neutral_Higgs_masses(Mh, dMh)
! Sets the theoretical mass uncertainty of the Higgs bosons.
!------------------------------------------------------------
 use usefulbits, only: theo,np,Hneut
 
 implicit none
 double precision,intent(out) :: Mh(np(Hneut)), dMh(np(Hneut))

 if(.not.allocated(theo))then
  stop 'No model information given!'
 endif
 if(np(Hneut).eq.0)then
  write(*,*)'Cannot access the neutral Higgs boson masses'
  write(*,*)'because np(Hneut) == 0.'
  stop 'error in subroutine get_neutral_Higgs_masses'
 endif

 Mh = theo(1)%particle(Hneut)%M 
 dMh = theo(1)%particle(Hneut)%dM
 
end subroutine get_neutral_Higgs_masses
!------------------------------------------------------------
subroutine finish_HiggsSignals
! This subroutine needs to be called right at the end, to close files
! and deallocate arrays
!------------------------------------------------------------
 use usefulbits, only : deallocate_usefulbits,debug,theo,debug,inputsub, &
   & file_id_debug1,file_id_debug2
 use S95tables, only : deallocate_Exptranges
 use theory_BRfunctions, only : deallocate_BRSM
 use datatables, only : deallocate_observables
 use usefulbits_HS, only : deallocate_usefulbits_HS, analyses
 use mc_chisq, only : deallocate_mc_observables
 use store_pathname_HS
 
!#if defined(NAGf90Fortran)
! use F90_UNIX_IO, only : flush
!#endif
      
 if(debug)then
  close(file_id_debug2)
  close(file_id_debug1)
 endif

 if(debug) write(*,*)'finishing off...'                      ; call flush(6)
 if(.not.allocated(theo))then
!  stop 'HiggsBounds_initialize  should be called first'
 if(debug) write(*,*) "HiggsBounds/HiggsSignals internal structure already deallocated!"
 else
  call deallocate_BRSM
  call deallocate_Exptranges
  call deallocate_usefulbits
  if (allocated(inputsub)) deallocate(inputsub)
 endif 
 call deallocate_mc_observables
 call deallocate_observables
 if(allocated(analyses)) deallocate(analyses)
 call deallocate_usefulbits_HS
! call system('rm -f '//trim(adjustl(pathname_HS))//'Expt_tables/analyses.txt')
 call system('rm -f HS_analyses.txt')
  
 if(debug) write(*,*)'finished'                              ; call flush(6)
 
end subroutine finish_HiggsSignals
!------------------------------------------------------------
! SOME HANDY WRAPPER SUBROUTINES
!------------------------------------------------------------
subroutine initialize_HiggsSignals_for_Fittino(nHiggsneut,nHiggsplus)
!------------------------------------------------------------
! Wrapper subroutine to intitialize HiggsSignals with the experimental
! dataset "latestresults", avoiding to specify this via a string argument.
!------------------------------------------------------------
 implicit none
 !--------------------------------------input
 integer,intent(in) :: nHiggsneut
 integer,intent(in) :: nHiggsplus
! character(LEN=19) :: Expt_string
 character(LEN=33) :: Expt_string 
! Expt_string = "Moriond2013_Fittino"
 Expt_string = "latestresults_April2013_inclusive"
 
 call initialize_HiggsSignals(nHiggsneut,nHiggsplus,Expt_string)
 
end subroutine initialize_HiggsSignals_for_Fittino
!------------------------------------------------------------
subroutine get_number_of_observables_wrapper(ntotal, npeakmu, npeakmh, nmpred, nanalyses)
!------------------------------------------------------------
 use io, only : get_number_of_observables
 
 implicit none
 integer, intent(out) :: ntotal, npeakmu, npeakmh, nmpred, nanalyses
 
 call get_number_of_observables(ntotal, npeakmu, npeakmh, nmpred, nanalyses) 
end subroutine get_number_of_observables_wrapper
!------------------------------------------------------------
subroutine get_ID_of_peakobservable_wrapper(ii, ID)
!------------------------------------------------------------
 use io, only : get_ID_of_peakobservable 
 implicit none
 integer, intent(in) :: ii
 integer, intent(out) :: ID

 call get_ID_of_peakobservable(ii, ID)
end subroutine get_ID_of_peakobservable_wrapper
!------------------------------------------------------------
subroutine get_peakinfo_from_HSresults_wrapper(obsID, mupred, domH, nHcomb)
!--------------------------------------------------------------------
 use io, only : get_peakinfo_from_HSresults
 
 implicit none
 integer, intent(in) :: obsID
 double precision, intent(out) :: mupred
 integer, intent(out) :: domH, nHcomb

 call get_peakinfo_from_HSresults(obsID, mupred, domH, nHcomb)
end subroutine get_peakinfo_from_HSresults_wrapper
!------------------------------------------------------------
subroutine print_cov_mh_to_file_wrapper(Hindex)
!------------------------------------------------------------
 use pc_chisq, only : print_cov_mh_to_file

 implicit none
 integer, intent(in) :: Hindex
 
 call print_cov_mh_to_file(Hindex)
end subroutine print_cov_mh_to_file_wrapper
!------------------------------------------------------------
subroutine print_cov_mu_to_file_wrapper
!------------------------------------------------------------
 use pc_chisq, only : print_cov_mu_to_file

 implicit none
 call print_cov_mu_to_file
end subroutine print_cov_mu_to_file_wrapper
!------------------------------------------------------------
