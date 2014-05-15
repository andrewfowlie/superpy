!******************************************************
! SuperPy, but it based on the
! example program HBandHSwithSLHA.f90
!******************************************************
program HBandHSwithSLHA
!
! In this example we run both HiggsBounds and HiggsSignals simultaneously
! on SLHA file(s). NOTE: The feature of selecting different experimental
! data in HiggsBounds is not fully supported here. In particular, the onlyL
! option of HiggsBounds does not work (it will be turned into LandH).
! If you want to take into account only LEP limits in HiggsBounds, you have
! to run both programs separately, using e.g. the example programs
! HBwithSLHA and HSwithSLHA.
!
! The SLHA file(s) has to contain the two HiggsBounds SLHA
! input blocks
!  HiggsBoundsInputHiggsCouplingsBosons
!  HiggsBoundsInputHiggsCouplingsFermions
! (see HiggsBounds (version 3 or more) manual for more details)
!
! Run with
!   ./HSwithSLHA npoints <stem>
! where npoints is the number of parameter points you would like to
! look at and each parameter point has a corresponding SLHA file
! e.g. the corresponding SLHA for the 5th point should be found at
!   <stem>.5
!
! Output
! The block HiggsSignalsResults will be added to each SLHA file.
!
!******************************************************

use io, only : HiggsSignals_SLHA_output, get_peakinfo_from_HSresults

implicit none

integer :: nH,nHplus,ndf,iter,HBresult,chan,ncombined
double precision :: obsratio
double precision :: Chisq, Chisq_mu, Chisq_mh, Pvalue
double precision, allocatable :: dMh(:), dCS(:), dBR(:)
integer :: i,npoints
integer,parameter :: fileid=78, fileid2=79
character(len=8) :: istring
character(len=300) :: inputfilename,outputfilename
character(len=300) :: stem
character(LEN=300) :: temp, tmpstring
integer :: number_args, ios

nH=3
nHplus=1

allocate(dMh(nH))

! Estimates of (relative!) systematic uncertainties:
!  dCS(1) - singleH       dBR(1) - gamma gamma
!  dCS(2) - VBF           dBR(2) - W W
!  dCS(3) - HW            dBR(3) - Z Z
!  dCS(4) - HZ            dBR(4) - tau tau
!  dCS(5) - ttH           dBR(5) - b bbar
allocate(dCS(5),dBR(5))

! The rate uncertainties for your model (typical MSSM estimates):
! SuperPy - unchanged from defaults for MSSM.
dCS = (/ 0.20D0, 0.028D0, 0.037D0, 0.051D0, 0.12D0 /)
dBR = (/ 0.054D0, 0.048D0, 0.048D0, 0.061D0, 0.028D0 /)

! Theoretical uncertainties on Higgs masses dMh (in GeV) for h,H,A:
! SuperPy - unchanged from defaults for MSSM.
dMh = (/ 2.0D0, 2.0D0, 0.0D0 /)

!------------------------------------------------------------------------------!
! Read arguments

! Read argument from comand line.
number_args = IARGC()
if(number_args .ne. 1) stop "Error: use ./SuperPy <SLHA FILE>"

! Read arguments into text strings.
call GETARG(1,inputfilename)
inputfilename = trim(inputfilename)
write(*,*) "SLHA file is:", inputfilename

!------------------------------------------------------------------------------!
! HiggsSignals options.

! Choose dataset -
! LandH = LEP, Tevatron and LHC data
call initialize_HiggsBounds(nH,nHplus,'LandH')

! Name of analysis folder.
call initialize_HiggsSignals(nH,nHplus,"latestresults")


! Output level (0: silent, 1: screen output, 2: even more output, ...
call setup_output_level(1)

! Higgs mass parametrization (1: box, 2:gaussian, 3:box+gaussian)
call setup_pdf(2)

! Number of iterations to find the (best) Higgs-to-peaks assignment
call setup_Higgs_to_peaks_assignment_iterations(0)

! If the mass-centered chi^2 method is used, can specify the dm_theory treatment
! 1: mass variation, 2: smearing of mu-plot with mass pdf
call setup_mcmethod_dm_theory(1)

!------------------------------------------------------------------------------!

! Read SLHA.
call HiggsBounds_input_SLHA(inputfilename)

!------------------------------------------------------------------------------!

! We want to use the mass variation treatment for the theoretical uncertainty.
! We set the mass uncertainties here (neutral Higgses mass errors
! set to dMh, charged Higgs mass error set to 0.0D0)

! SuperPy - NB careful not to clash with call setup_mcmethod_dm_theory(1)
! I don't understand the program structure.
call HiggsBounds_set_mass_uncertainties(dMh,0.0D0)

!------------------------------------------------------------------------------!

! Run HiggsBounds
call run_HiggsBounds(HBresult, chan, obsratio, ncombined)

!------------------------------------------------------------------------------!

! Read SLHA.
! SuperPy - again?! Perhaps variables are overwritten/discarded in HiggsBounds.
call HiggsBounds_input_SLHA(inputfilename)

!------------------------------------------------------------------------------!

! Pass the Higgs mass uncertainty to HiggsSignals
call HiggsSignals_neutral_input_MassUncertainty(dMh)
! Set the production and decay rate uncertainties for the model
call setup_rate_uncertainties(dCS, dBR)

!------------------------------------------------------------------------------!
! Run HiggsSignals
call run_HiggsSignals(1, Chisq_mu, Chisq_mh, Chisq, ndf, Pvalue)

!------------------------------------------------------------------------------!
! Tidy up.
call finish_HiggsBounds
call finish_HiggsSignals
end program