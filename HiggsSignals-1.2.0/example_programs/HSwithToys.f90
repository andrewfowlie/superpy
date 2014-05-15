!--------------------------------------------------------------------------------------
program HSwithToys
!
! This example program is part of HiggsSignals (TS 25/01/2013).
!--------------------------------------------------------------------------------------
! This example program shows how HiggsSignals (HS) can be run on toy experiments. We
! first run HS with SM input with a Higgs mass of 126 GeV. From the result of this
! run we can obtain the predicted signal strength modifiers for each peak observable
! directly from HiggsSignals. As a demonstration, these are then set as observed (toy)
! signal rates of the various observables. We then re-run HiggsSignals on the new
! observables on the same model, i.e. observed and predicted signal rates are equal,
! thus the resulting chi^2 should be zero.
!--------------------------------------------------------------------------------------
 use theory_colliderSfunctions
 use usefulbits, only : vsmall 
 use usefulbits_hs,only : analyses
 use pc_chisq, only : print_peaks_to_LaTeX
 use io, only : get_peakinfo_from_HSresults, get_number_of_observables,&
 &              get_ID_of_peakobservable, HiggsSignals_create_SLHA_output,&
 &              HiggsSignals_create_SLHA_output_default
  	  
 implicit none
 integer :: nH,nHplus,ndf, ii, jj, kk
 double precision :: Chisq, Chisq_mu, Chisq_mh, Pvalue
 double precision, allocatable :: dMh(:)
 double precision :: dCS(5),dBR(5),dggh, dbbh
  double precision :: SMGamma_h
  double precision :: Mh,GammaTotal,g2hjss_s,g2hjss_p,g2hjcc_s,g2hjcc_p,&
&                     g2hjbb_s,g2hjbb_p,g2hjtt_s,g2hjtt_p,&
&                     g2hjmumu_s,g2hjmumu_p,g2hjtautau_s,g2hjtautau_p,&
&                     g2hjWW,g2hjZZ,g2hjZga,g2hjgaga,g2hjgg,g2hjggZ,&
&                     g2hjhiZ,BR_hjhihi,BR_hjinvisible  
  double precision :: mupred
  integer :: domH, nHcomb
  integer :: ntotal, npeakmu, npeakmh, nmpred, nanalyses, ID
  integer :: i,npoints
  character(len=8) :: istring
  character(len=300) :: inputfilename,outputfilename
  character(len=300) :: stem
  character(LEN=300) :: temp
  integer :: number_args, stat

!---Note here: We only run HiggsSignals on the lightest Higgs boson. This can be easily
!---extended to all 3 MSSM neutral Higgs bosons. In that case, the effective couplings
!---and mass uncertainties have to be given as arrays of size=nH (Cf. the Higgsbounds
!---manual for HB-3.x.x for how to call HiggsBounds_neutral_input_effC correctly!)
  nH=1
  nHplus=0

  allocate(dMh(nH))
!--Set the (relative!) rate uncertainties for your model:
!  dCS(1) - singleH				dBR(1) - gamma gamma
!  dCS(2) - VBF					dBR(2) - W W
!  dCS(3) - HW					dBR(3) - Z Z
!  dCS(4) - HZ					dBR(4) - tau tau
!  dCS(5) - ttH					dBR(5) - b bbar
!--Here, we first set them to the SM. Later, we determine the uncertainty of singleH
!--production, dCS(1), from the effective couplings. SM values are taken from
!--https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CrossSections
! dCS = (/ 0.147D0, 0.028D0, 0.037D0, 0.051D0, 0.12D0 /)
! dBR = (/ 0.054D0, 0.048D0, 0.048D0, 0.061D0, 0.028D0 /)
!!--Set the rate uncertainties of gluon-gluon fusion and bb->H here:
! dggh = 0.147D0
! dbbh = 0.200D0
!--n.b. have to set theoretical uncertainties on Higgs mass dMh (in GeV):
  dMh = (/ 0.0D0 /)
!-------------------------- HiggsSignals ------------------------------!

!---- Initialize HiggsSignals and pass the name of the experimental analysis folder  ----!
  call initialize_HiggsSignals(nH,nHplus,"latestresults") 
!  call initialize_HiggsSignals(nH,nHplus,"latestresults")   
!  call setup_thu_observables(1)   
!---- Set the Higgs mass parametrization (1: box, 2:gaussian, 3:box+gaussian)	 	 ----!
  call setup_pdf(2)
!---- Set the output level (0: silent, 1: screen output, 2: even more output,...) 	 ----!
  call setup_output_level(0)
!---- Pass the Higgs mass uncertainty to HiggsSignals							 	 ----!
  call HiggsSignals_neutral_input_MassUncertainty(dMh)
!---- Use symmetric rate errors? (0: original(default), 1: averaged-symmetrical)     ----!
! call setup_symmetricerrors(0)
!---- Allow anti-correlated signal strength measurements? (0: no, 1: yes(default) )  ----!
! call setup_anticorrelations_in_mu(1)
!---- Setup a wider assignment range                                                 ----!
 call setup_assignmentrange(2.0D0)
! call setup_assignmentrange_massobservables(2.0D0)

!----HiggsBounds/Signals effective couplings input. 
!  	 These have to be inserted for the model which we want to test, i.e. we would have
!    to write an interface to set via arguments in the executables call, or reading
!    in a text file, etc.
!----For now, we set them by hand to the SM values (for demonstration):
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
  Mh=dble(125.8)
  GammaTotal=SMGamma_h(Mh)
	  
!-Calculate theoretical uncertainties of singleH production from ggh and bbh effective couplings. 
!  dCS(1) = get_singleH_uncertainty(dggh, dbbh, g2hjgg, g2hjbb_s+g2hjbb_p, mh)	
!  call setup_rate_uncertainties(dCS, dBR)

!-Set the HiggsSignals input
  call HiggsBounds_neutral_input_effC(Mh,GammaTotal,&
&    g2hjss_s,g2hjss_p,g2hjcc_s,g2hjcc_p,g2hjbb_s,g2hjbb_p,&
&    g2hjtt_s,g2hjtt_p,&
&    g2hjmumu_s,g2hjmumu_p,g2hjtautau_s,g2hjtautau_p,&
&    g2hjWW,g2hjZZ,g2hjZga,g2hjgaga,g2hjgg,g2hjggZ,&
&    g2hjhiZ, BR_hjinvisible,BR_hjhihi)

!-Run HS on the original experimental data in order to evaluate the model predictions	  
  call run_HiggsSignals(1, Chisq_mu, Chisq_mh, Chisq, ndf, Pvalue)
!-Print out the observables to a LaTeX table
!  call print_peaks_to_LaTeX

!-Get the number of the peak-observables (Don't care about ntotal, npeakmh, nmpred, nanalyses)
  call get_number_of_observables(ntotal, npeakmu, npeakmh, nmpred, nanalyses)

!-We now want to set the measurements to those values predicted by the model.
!-The mass measurement for each peak observable will be set to Mh here.

!-Loop over the number of peak observables
  do kk=1,npeakmu
!--Get for each peak observable its unique ID:
   call get_ID_of_peakobservable(kk, ID)
!--Get the predicted signal strength modifier (mupred) for this peak observable:
   call get_peakinfo_from_HSresults(ID, mupred, domH, nHcomb)
!--Assign this value as (toy) measurement for this peak observable:
   call assign_toyvalues_to_peak(ID, mupred, Mh)
  enddo
  
  call setup_output_level(0)	!-Do a print-out to the screen
    
!-Set the HiggsSignals input (again!)
  call HiggsBounds_neutral_input_effC(Mh,GammaTotal,&
&    g2hjss_s,g2hjss_p,g2hjcc_s,g2hjcc_p,g2hjbb_s,g2hjbb_p,&
&    g2hjtt_s,g2hjtt_p,&
&    g2hjmumu_s,g2hjmumu_p,g2hjtautau_s,g2hjtautau_p,&
&    g2hjWW,g2hjZZ,g2hjZga,g2hjgaga,g2hjgg,g2hjggZ,&
&    g2hjhiZ, BR_hjinvisible,BR_hjhihi)

!-Now, we run on the toy observables with the same input, i.e. model predictions and
!-measurements are equal and thus the chi^2 should be zero.	  
  call run_HiggsSignals(1, Chisq_mu, Chisq_mh, Chisq, ndf, Pvalue)

!-Create a new SLHA file with the HiggsSignals output blocks.
!-The second argument controls how much is written
!  (0: only the BLOCK 'HiggsSignalsResults', 1: full HiggsSignals SLHA output)
!-The new file must not exist:
! (note: this system call does not work with ifort)
!  call system('rm -f results/HSwithToys.slha',status=stat)
  call HiggsSignals_create_SLHA_output("results/HSwithToys.slha",0)
!-Alternatively, we could use
!  call HiggsSignals_create_SLHA_output_default(0)
!-where the filename is set to "HS-output.slha".

  call finish_HiggsSignals
 
 contains
!--------------------------------------------------------------------------------------
 function get_singleH_uncertainty(dggh, dbbh, g2hgg, g2hbb, mh)
 ! This function evaluates the uncertainty of single Higgs production from an interpolation
 ! of the uncertainties on gg->H and bb->H, using the squared effective couplings g2hgg and
 ! g2hbb. It uses internal HiggsBounds functions from the module theory_colliderSfunctions.
!--------------------------------------------------------------------------------------

 double precision, intent(in) :: dggh, dbbh, g2hgg, g2hbb, mh
 double precision :: get_singleH_uncertainty
 
 if(g2hgg.le.vsmall.and.g2hbb.le.vsmall) then
  get_singleH_uncertainty = 0.0D0
 else 
  get_singleH_uncertainty = ( g2hgg*LHC8_rH_gg(mh)*dggh + g2hbb*LHC8_rH_bb(mh)*dbbh )/ &
 & 	                        ( g2hgg*LHC8_rH_gg(mh)      + g2hbb*LHC8_rH_bb(mh)      )
 endif
 
 end function get_singleH_uncertainty
!--------------------------------------------------------------------------------------
end program HSwithToys
!--------------------------------------------------------------------------------------
