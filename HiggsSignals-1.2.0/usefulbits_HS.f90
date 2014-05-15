!--------------------------------------------------------------------
! This file is part of HiggsSignals (TS 03/03/2013)
!--------------------------------------------------------------------
module usefulbits_HS
!--------------------------------------------------------------------
 implicit none

 character(LEN=9),parameter :: HSvers='1.2.0'
 integer,parameter :: f_dmh=94
 character(LEN=4) :: runmode
 character(LEN=100) :: Exptdir
!--------------------------------------------------------------------
!--IN TESTING PHASE: (do not change this in official version)
 logical :: withcorrexpsyst = .True. !(correlated experimental systematics)
 logical :: additional_output = .False.
!--END 
!------------------------Control parameters -------------------------
! Values can be changes with specific user subroutines.
 logical :: usetoys = .False.
 logical :: usescalefactor = .False.
 logical :: useSMweights = .False.
 logical :: correlations_mu = .True.
 logical :: correlations_mh = .True. 
 logical :: minimalchisq = .False.
 logical :: maximalchisq = .False.
 logical :: useSMtest = .False.
 logical :: SLHAdetailed = .True.
 logical :: newSLHAfile = .False.
 logical :: symmetricerrors = .False.
 logical :: anticorrmu = .True.
 logical :: anticorrmh = .True.
 logical :: absolute_errors = .False.
 logical :: THU_included = .True.

! CAREFUL: Only switch this on if you know what you are doing. This is still in
!          testing phase!
 logical :: normalize_rates_to_reference_position = .False.

 double precision :: eps 
 double precision :: assignmentrange = 1.0D0  ! This gives the mass range
                                              ! (in standard deviations), in which
                                              ! the Higgs is forced to be assigned to
                                              ! a peak observable.
 double precision :: assignmentrange_massobs = 1.0D0  
                                              ! This gives the mass range
                                              ! (in standard deviations), in which
                                              ! the Higgs is forced to be assigned to
                                              ! peak observables, which have a mass
                                              ! measurement.

 integer :: output_level = 0
 integer :: iterations = 0		! default value: 0
 								! 1 -> try to assign as many Higgs bosons as possible to
 								! the observable, Higgs-to-peak assignment is based on
 								! Higg mass covariance matrices with maximal
 								! correlations.
 								! >1 -> use the covariance matrix of previous iteration.
 integer :: pdf	= 2				! default value: 2
 								! will automatically be set to 2 if not changed by the user
 								! via using subroutine set_pdf before.
 								! (1,2,3) = (box, gaussian, theory-box + exp-gaussian)
 integer :: Nparam = 0			! Number of free model parameters (entering the Pvalue)
 								! Can be specified directly here or via the subroutine
 								! setup_nparam 								
!-------------------------------------------------------------------- 								
 integer :: nanalys 			!Total number of relevant analyses                       
 double precision, parameter :: vlarge=1000000000000.0D0
!--------------------- Default rate uncertainties -------------------

 type rate_uncertainties 
!- dCS_SM and dBR_SM for the SM 
!- (from LHC HXSWG Yellow Report 3, arXiv:1307.1347)
!- dCS and dBR hold the model's rate uncertainties. Can be changed by user
!- with subroutine setup_rate_uncertainties. Default values are those of the SM.
  double precision :: dCS_SM(5) = (/ 0.147D0, 0.028D0, 0.037D0, 0.060D0, 0.12D0 /)  
  double precision :: dCS(5) = (/ 0.147D0, 0.028D0, 0.037D0, 0.060D0, 0.12D0 /)  
!  double precision :: dBR_SM(5) = (/ 0.054D0, 0.048D0, 0.048D0, 0.061D0, 0.028D0 /)
!  double precision :: dBR(5) = (/ 0.054D0, 0.048D0, 0.048D0, 0.061D0, 0.028D0 /)
!- EDIT (TS 21/06/2013): Add new decay modes:
!- Channels: gammagamma, WW, ZZ, tautau, bb, Zgamma, cc, mumu, gg
  double precision :: dBR_SM(9) = (/ 0.054D0, 0.048D0, 0.048D0, 0.061D0, 0.028D0,&
&                                    0.090D0, 0.122D0, 0.060D0, 0.100D0 /)  
  double precision :: dBR(9) = (/ 0.054D0, 0.048D0, 0.048D0, 0.061D0, 0.028D0,&
&                                    0.090D0, 0.122D0, 0.060D0, 0.100D0 /)
!--- IMPORTANT NOTE:
!-
!- The arrays dCS_SM, dCS, dBR_SM, dBR have been introduced in HiggsSignals-1.0.0
!- to hold the estimated theoretical uncertainties. These do not include correlations
!- via parametric uncertainties (e.g. scale, PDFs,...) or correlations in the BRs introduced
!- by the uncertainty of the total widths.
!-
!- Since HiggsSignals-1.1.0 the theoretical uncertainties for the cross sections and
!- branching ratios are evaluated with toy MC scripts including the correlations of
!- parametric error sources. The resulting covariance matrices are included per default
!- if the files "BRcov.in" and "XScov.in" are present in the main HiggsSignals directory.
!- If not, HiggsSignals will give a warning and use the old method.
!- The covariance matrices can also be re-evaluated by the user with the scripts
!- "smearErrorsBR.cpp" and "smearErrorsXS.cpp", which can be found in the directory
!- <HiggsSignals-main-directory>/supplements/
!-
!---
  logical :: BRcov_ok=.False.
  logical :: CScov_ok=.False.
  logical :: usecov =.True.
  
  double precision, dimension(9,9) :: BRcovSM = 0.0D0
  double precision, dimension(9,9) :: BRcov = 0.0D0
  double precision, dimension(5,5) :: CScovSM = 0.0D0
  double precision, dimension(5,5) :: CScov = 0.0D0
    
!--- ILC cross section uncertainties (under development)  
!--- (none, none, WBF, ZH, ttH)
  double precision :: dCS_ILC_SM(5) = (/ 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0 /)  
  double precision :: dCS_ILC(5) = (/ 0.0D0, 0.0D0, 0.01D0, 0.005D0, 0.01D0 /)  


 end type  
 type(rate_uncertainties), save :: delta_rate 
 
!-------------- Type definitions of internal structures --------------
 type neutHiggs
  double precision :: m, dm, mu
  integer :: mp_test	! This variable is set to 1 (0) if the Higgs is (not) being tested in the m-pred chi^2 method.
  integer :: id
 end type
 !-Will contain info about all neutral Higgs for every considered observable, i.e.
!-neutHiggses has dimensions (number(observables),nH)
 type(neutHiggs), allocatable :: neutHiggses(:,:)

type mutable
  integer :: id,nx,particle_x !see usefulbits.f90 for key to particle codes n.b. they're NOT pdg
  character(LEN=45) :: label
  character(LEN=100) :: desc
  character(LEN=3) :: expt  
  character(LEN=10) :: collider 
  character(LEN=10) :: collaboration   
  double precision :: lumi,dlumi,energy			! dlumi in %
!--TESTING correlated experimental systematics:
!  double precision, dimension(4) :: correxpsyst
!--END  
  double precision :: xmax,xmin,sep,deltax
  double precision :: deltam
  character(LEN=100) :: assignmentgroup
  integer :: mhchisq
  double precision, allocatable :: mass(:)
  double precision, allocatable :: mu(:,:) ! in mu(a,b), a=row, b=1,2,3 for low,obs,up 
  integer :: Nc								  ! Number of channels
  integer, allocatable :: channel_id(:)		  ! Considered channels array, dim(Nc)
  character(LEN=10),allocatable :: channel_description(:,:)  
  double precision, allocatable :: channel_eff(:) 	  ! Channel efficiencies, dim(Nc) 
  double precision, allocatable :: channel_eff_ratios(:) 	  ! Channel efficiency ratios (model vs. SM), dim(Nc)   
  double precision, allocatable :: channel_w(:,:) ! Channel weights, dim(Nc, NHiggs)
  double precision, allocatable :: channel_w_corrected_eff(:,:) ! Channel weights, dim(Nc, NHiggs)  
  double precision, allocatable :: channel_systSM(:,:) ! Channel systematics of SM, dim(Nc, NHiggs)  
  double precision, allocatable :: channel_syst(:,:) ! Channel systematics, dim(Nc, NHiggs)
  double precision, allocatable :: channel_mu(:,:) ! SM normalized channel rates, dim(Nc, NHiggs)											  
  double precision :: eff_ref_mass			  ! Reference Higgs mass for quoted efficiency
  integer :: npeaks
  double precision, allocatable :: Toys_mhobs(:)
  double precision, allocatable :: Toys_muobs(:)
  double precision :: scale_mu
 end type

 type mupeak
  integer :: id      
  integer :: ilow,iup,ipeak
  double precision :: mpeak
  double precision :: dm
  double precision :: mu
  double precision :: mu_original  
  double precision :: scale_mu!, scale_mh  
  double precision :: dmuup,dmulow			  ! Upper and lower cyan band 	
  double precision :: dmuup0sq, dmulow0sq	  ! Cyan band squared subtracted by correlated uncertainties
  !-Peak object should contain everything needed for constructing the covariance matrices
  integer :: Nc								  ! Number of channels
  integer, allocatable :: channel_id(:)		  ! Considered channels array, dim(Nc)
  double precision, allocatable :: channel_eff(:) 	  ! Channel efficiencies, dim(Nc)   
  integer, allocatable :: Higgs_comb(:)		  ! Assigned Higgs combination, dim(NHiggs)
  character(LEN=100) :: assignmentgroup  
  type(neutHiggs), allocatable :: Higgses(:)
  integer :: domH							  !	index of dominantly contributing Higgs
  integer :: NHiggs_comb					  ! Number of combined Higgses
  integer :: Higgs_assignment_forced
  integer :: undo_assignment
  !--These arrays contain only the information about all Higgs bosons
  !--(need to have this for every peak separately because it can depend on the efficiencies
  !-- which are given for each peak separately)
  double precision, allocatable :: channel_w_allH(:,:)      ! Channel weights, dim(Nc, NHiggs)
  double precision, allocatable :: channel_w_corrected_eff_allH(:,:) ! Channel weights with corrected efficiencies, dim(Nc, NHiggs)  
  double precision, allocatable :: channel_systSM_allH(:,:) ! Channel systematics of SM, dim(Nc, NHiggs)  
  double precision, allocatable :: channel_syst_allH(:,:)   ! Channel systematics, dim(Nc, NHiggs)
  double precision, allocatable :: channel_mu_allH(:,:)     ! SM normalized channel rates, dim(Nc, NHiggs)											  	
  !--These arrays contain only the information about the chosen Higgs combination:
  double precision, allocatable :: channel_w(:)             ! Channel weights, dim(Nc)
  double precision, allocatable :: channel_w_corrected_eff(:) ! Channel weights with corrected efficencies, dim(Nc)
  double precision, allocatable :: channel_systSM(:)        ! Channel systematics, dim(Nc)  
  double precision, allocatable :: channel_syst(:)          ! Channel systematics, dim(Nc)
  double precision, allocatable :: channel_mu(:)            ! SM normalized channel rates, dim(Nc)
  double precision, allocatable :: channel_w_model(:)
  double precision :: total_mu
  double precision :: dlumi
  !-- Chisq values (mu and mh parts, total) after taking into account correlations with
  !-- other peaks:
  double precision :: chisq_mu
  double precision :: chisq_mh  
  double precision :: chisq_tot
  double precision :: chisq_max
  integer :: internalnumber
 end type

type mp_neutHiggs
!-This object is a Higgs or Higgscluster which are separately
!-tested with the predicted mass chi^2 method.
 type(neutHiggs), allocatable :: Higgses(:)
 double precision :: m, dm, mu
 integer :: mp_test
 double precision :: mu_obs, dmu_low_obs, dmu_up_obs, dmu_low0_obs, dmu_up0_obs, m_obs
 double precision, allocatable :: channel_w_model(:)
 double precision, allocatable :: channel_mu(:)
 double precision :: total_mu
 integer :: Higgscomb
 integer :: domH
 double precision :: chisq 
!-n.b. these are the smeared observed signal strengths for this Higgs boson
end type

type mpred
 type(mp_neutHiggs), allocatable :: mp_Higgses(:) 
 double precision :: mupred
end type

type observable
 integer :: id
 integer :: obstype
 type(mupeak) :: peak
 type(mutable) :: table
 type(neutHiggs), allocatable :: Higgses(:)
end type
type(observable), allocatable :: obs(:)

type tablelist
 integer :: Npeaks
 integer :: id
 type(mutable) :: table
 type(mupeak), allocatable :: peaks(:)
 type(mpred) :: mpred
 type(neutHiggs), allocatable :: Higgses(:)	
	! This object holds primarily the Higgs boson predictions
 	! corresponding to tablelist%table. It corresponds to the full 
 	! muplot if it is implemented (to enable the mpred-method).
end type
type(tablelist), allocatable :: analyses(:)

type HSresults                  
 double precision :: Pvalue = -1.0D0           
 double precision :: Chisq
 double precision :: Chisq_mu, Chisq_mpred, Chisq_peak_mu
 double precision :: Chisq_mh
 double precision, allocatable :: mupred(:)
 integer, allocatable :: domH(:)    
 integer, allocatable :: nH(:)  
 integer, allocatable :: obsID(:)  
 integer :: nobs, nobs_mpred, nobs_peak_mu, nobs_peak_mh, nanalysis
end type            
type(HSresults), allocatable :: HSres(:)

!----------------------- Covariance matrices ----------------------
 double precision, allocatable :: cov(:,:)
 double precision, allocatable :: cov_mhneut(:,:,:)
 double precision, allocatable :: cov_mhneut_max(:,:,:)
 double precision, allocatable :: cov_mp(:,:)
 double precision, allocatable :: cov_mu_tot(:,:)
 double precision, allocatable :: mu_vector(:)
!--------------------------------------------------------------------  
 contains
!--------------------------------------------------------------------
 subroutine HiggsSignals_info
!--------------------------------------------------------------------
  implicit none

  write(*,*)
  write(*,*)"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
  write(*,*)"~                                                        ~"
  write(*,*)"~      	      HiggsSignals "//adjustl(HSvers)//"	          ~"
  write(*,*)"~                                                        ~"
  write(*,*)"~      Philip Bechtle, Sven Heinemeyer, Oscar Stål,      ~"
  write(*,*)"~              Tim Stefaniak, Georg Weiglein             ~"
  write(*,*)"~                                                        ~"
  write(*,*)"~                    arXiv:1305.1933                     ~"
  write(*,*)"~                                                        ~"
  write(*,*)"~ It is based on the HiggsBounds-4 Fortran library.      ~"
  write(*,*)"~ Please consult and cite also the following references  ~"
  write(*,*)"~ for the HiggsBounds program                            ~"
  write(*,*)"~                                                        ~"  
  write(*,*)"~    arXiv:0811.4169, arXiv:1102.1898, arXiv:1311.0055   ~"  
  write(*,*)"~                                                        ~"    
  write(*,*)"~            http://higgsbounds.hepforge.org             ~"
  write(*,*)"~                                                        ~"
  write(*,*)"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
  write(*,*)
  write(*,*)" HiggsSignals collects together results from "
  write(*,*)
  write(*,*)"    * the ATLAS and CMS Collaborations"
  write(*,*)"    * the CDF and D0 Collaborations"
  write(*,*)"    * the program HDECAY (arXiv:hep-ph/9704448)"  
  write(*,*)"    * LHC Higgs Cross Section Working Group"
  write(*,*)"      (arXiv:1101.0593, arXiv:1201.3084, "
  write(*,*)"       arXiv:1307.1347 and ref. therein)"
  write(*,*)
  
  
 end subroutine HiggsSignals_info
!--------------------------------------------------------------------
 subroutine print_dble_matrix(mat, title)
 !--------------------------------------------------------------------
  implicit none
  double precision, dimension(:,:), intent(in) :: mat(:,:)
  character(LEN=50), intent(in), optional :: title
  integer :: i
  
  if(present(title)) then
   write(*,*)"#*************************************************************************#"
   write(*,*)"# ",trim(title)
  endif 
  write(*,*) "#*************************************************************************#"
  do i=lbound(mat,dim=1),ubound(mat,dim=1)
	write(*,*) mat(i,:)
  enddo	
  write(*,*) "#*************************************************************************#"
 end subroutine print_dble_matrix
!--------------------------------------------------------------------
 subroutine deallocate_usefulbits_HS
!--------------------------------------------------------------------  
  implicit none
  integer :: i
!  deallocate(neutHiggses)
  if(allocated(HSres)) then
   do i=lbound(HSres, dim=1), ubound(HSres, dim=1)
    if(allocated(HSres(i)%mupred)) deallocate(HSres(i)%mupred)
    if(allocated(HSres(i)%domH)) deallocate(HSres(i)%domH)
    if(allocated(HSres(i)%nH)) deallocate(HSres(i)%nH)   
   enddo 
   deallocate(HSres)
  endif	

  call deallocate_covariance_matrices
  
 end subroutine deallocate_usefulbits_HS
!--------------------------------------------------------------------  
 subroutine deallocate_covariance_matrices
!--------------------------------------------------------------------  
  implicit none
  
  if(allocated(cov)) deallocate(cov)
  if(allocated(cov_mhneut)) deallocate(cov_mhneut)
  if(allocated(cov_mhneut_max)) deallocate(cov_mhneut_max)
  if(allocated(cov_mp)) deallocate(cov_mp)
  if(allocated(cov_mu_tot)) deallocate(cov_mu_tot)
  if(allocated(mu_vector)) deallocate(mu_vector)
  
 end subroutine deallocate_covariance_matrices
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       StrCompress
!
! PURPOSE:
!       Subroutine to return a copy of an input string with all whitespace
!       (spaces and tabs) removed.
!
! CALLING SEQUENCE:
!       Result = StrCompress( String,  &  ! Input
!                             n = n    )  ! Optional Output
!
! INPUT ARGUMENTS:
!       String:         Character string to be compressed.
!                       UNITS:      N/A
!                       TYPE:       CHARACTER(*)
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
! OPTIONAL OUTPUT ARGUMENTS:
!       n:              Number of useful characters in output string
!                       after compression. From character n+1 -> LEN(Input_String)
!                       the output is padded with blanks.
!                       UNITS:      N/A
!                       TYPE:       INTEGER
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(OUT), OPTIONAL
!
! FUNCTION RESULT:
!       Result:         Input string with all whitespace removed before the
!                       first non-whitespace character, and from in-between 
!                       non-whitespace characters.
!                       UNITS:      N/A
!                       TYPE:       CHARACTER(LEN(String))
!                       DIMENSION:  Scalar
!
! EXAMPLE:
!       Input_String = '  This is a string with spaces in it.'
!       Output_String = StrCompress( Input_String, n=n )
!       WRITE( *, '( a )' ) '>',Output_String( 1:n ),'<'
!   >Thisisastringwithspacesinit.<
!
!       or
!
!       WRITE( *, '( a )' ) '>',TRIM( Output_String ),'<'
!   >Thisisastringwithspacesinit.<
!
! PROCEDURE:
!       Definitions of a space and a tab character are made for the
!       ASCII collating sequence. Each single character of the input
!       string is checked against these definitions using the IACHAR()
!       intrinsic. If the input string character DOES NOT correspond 
!       to a space or tab, it is not copied to the output string.
!
!       Note that for input that ONLY has spaces or tabs BEFORE the first
!       useful character, the output of this function is the same as the
!       ADJUSTL() instrinsic.
!
! CREATION HISTORY:
!       Written by:     Paul van Delst, CIMSS/SSEC 18-Oct-1999
!                       paul.vandelst@ssec.wisc.edu
!
!:sdoc-:
!--------------------------------------------------------------------
 FUNCTION StrCompress( Input_String, n ) RESULT( Output_String )
    ! Arguments
    CHARACTER(*),      INTENT(IN)  :: Input_String
    INTEGER, OPTIONAL, INTENT(OUT) :: n
    ! Function result
    CHARACTER(LEN(Input_String)) :: Output_String
    ! Local parameters
    INTEGER, PARAMETER :: IACHAR_SPACE = 32
    INTEGER, PARAMETER :: IACHAR_TAB   = 9
    ! Local variables
    INTEGER :: i, j
    INTEGER :: IACHAR_Character

    ! Setup
    ! -----
    ! Initialise output string
    Output_String = ' '
    ! Initialise output string "useful" length counter
    j = 0

    ! Loop over string contents character by character
    ! ------------------------------------------------
    DO i = 1, LEN(Input_String)

      ! Convert the current character to its position
      ! in the ASCII collating sequence
      IACHAR_Character = IACHAR(Input_String(i:i))

      ! If the character is NOT a space ' ' or a tab '->|'
      ! copy it to the output string.
      IF ( IACHAR_Character /= IACHAR_SPACE .AND. &
           IACHAR_Character /= IACHAR_TAB         ) THEN
        j = j + 1
        Output_String(j:j) = Input_String(i:i)
      END IF

    END DO

    ! Save the non-whitespace count
    ! -----------------------------
    IF ( PRESENT(n) ) n = j

  END FUNCTION StrCompress
!--------------------------------------------------------------------            
end module usefulbits_HS                        
!--------------------------------------------------------------------