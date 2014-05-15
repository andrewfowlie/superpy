!--------------------------------------------------------------------------------------
program HS_efficiencies
!
! This example program is part of HiggsSignals (TS 27/08/2013).
!--------------------------------------------------------------------------------------
! This example program demonstrates how models with different channel efficiencies 
! of the Higgs searches can be treated in HiggsSignals (HS).
!--------------------------------------------------------------------------------------
! use theory_colliderSfunctions
! use usefulbits, only : vsmall 
! use usefulbits_hs,only : analyses
! use pc_chisq, only : print_peaks_to_LaTeX
 use io, only : get_peakinfo_from_HSresults, get_number_of_observables,&
 &              get_ID_of_peakobservable, get_peak_channels
  	  
 implicit none
 integer :: nH,nHplus,ndf, jj, kk, ll, mm
 double precision :: Chisq, Chisq_mu, Chisq_mh, Pvalue
 double precision, allocatable :: dMh(:)
 double precision :: dCS(5),dBR(5),dggh, dbbh
  double precision :: SMGamma_h
  double precision :: Mh,GammaTotal,g2hjss_s,g2hjss_p,g2hjcc_s,g2hjcc_p,&
&                     g2hjbb_s,g2hjbb_p,g2hjtt_s,g2hjtt_p,&
&                     g2hjmumu_s,g2hjmumu_p,g2hjtautau_s,g2hjtautau_p,&
&                     g2hjWW,g2hjZZ,g2hjZga,g2hjgaga,g2hjgg,g2hjggZ,&
&                     g2hjhiZ,BR_hjhihi,BR_hjinvisible  
  double precision :: mupred, etaZ, sf, mu_average
  integer :: domH, nHcomb
  integer :: ntotal, npeakmu, npeakmh, nmpred, nanalyses, ID, pID
  integer :: i,npoints
  character(len=8) :: istring
  character(len=300) :: inputfilename,outputfilename
  character(len=300) :: stem
  character(LEN=300) :: temp
  integer :: number_args, stat
  integer :: Nc
  integer, allocatable :: channel_IDs(:)
  double precision, allocatable :: efficiencies(:),eff_ratio(:)

!---Note here: We only run HiggsSignals on the lightest Higgs boson. This can be easily
!---extended to all 3 MSSM neutral Higgs bosons. In that case, the effective couplings
!---and mass uncertainties have to be given as arrays of size=nH (Cf. the Higgsbounds
!---manual for HB-3.x.x for how to call HiggsBounds_neutral_input_effC correctly!)
  nH=1
  nHplus=0

  allocate(dMh(nH))
!--n.b. have to set theoretical uncertainties on Higgs mass dMh (in GeV):
  dMh = (/ 0.0D0 /)
!-------------------------- HiggsSignals ------------------------------!

!---- Initialize HiggsSignals and pass the name of the experimental analysis folder  ----!
  call initialize_HiggsSignals(nH,nHplus,"latestresults") 
!---- Set the Higgs mass parametrization (1: box, 2:gaussian, 3:box+gaussian)	 	 ----!
  call setup_pdf(2)
!---- Set the output level (0: silent, 1: screen output, 2: even more output,...) 	 ----!
  call setup_output_level(0)
!---- Pass the Higgs mass uncertainty to HiggsSignals							 	 ----!
  call HiggsSignals_neutral_input_MassUncertainty(dMh)
!---- Use symmetric rate errors? (0: original(default), 1: averaged-symmetrical)     ----!
 call setup_symmetricerrors(0)
!---- Allow anti-correlated signal strength measurements? (0: no, 1: yes(default) )  ----!
 call setup_anticorrelations_in_mu(0)
!---- Setup a wider assignment range                                                 ----!
 call setup_assignmentrange(10.0D0)

!----HiggsBounds/Signals effective couplings input. 
!  	 These have to be inserted for the model which we want to test, i.e. we would have
!    to write an interface to set via arguments in the executables call, or reading
!    in a text file, etc.
!----For now, we set them by hand to the SM values except for the vector boson couplings,
!    where we set them to a scale factor sf (for demonstration):

 do mm=1,4
  select case(mm)
   case(1)
    outputfilename="results/HS_efficiencies_kv0p5.dat"
   case(2)
    outputfilename="results/HS_efficiencies_kv1p0.dat"
   case(3)
    outputfilename="results/HS_efficiencies_kv1p5.dat"
   case(4)
    outputfilename="results/HS_efficiencies_kv2p0.dat"
   case default
   end select
  open(21,file=outputfilename)

  sf = 0.0D0+(mm)*0.5D0
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
  g2hjWW=sf
  g2hjZZ=sf
  g2hjZga=sf
  g2hjgaga=1d0
  g2hjgg=1d0
  g2hjggZ=sf
  g2hjhiZ=0d0
  BR_hjhihi=0d0
  BR_hjinvisible=0d0
  Mh=dble(125.8)
  GammaTotal=SMGamma_h(Mh)
	  
!-Set the HiggsSignals input
  call HiggsBounds_neutral_input_effC(Mh,GammaTotal,&
&    g2hjss_s,g2hjss_p,g2hjcc_s,g2hjcc_p,g2hjbb_s,g2hjbb_p,&
&    g2hjtt_s,g2hjtt_p,&
&    g2hjmumu_s,g2hjmumu_p,g2hjtautau_s,g2hjtautau_p,&
&    g2hjWW,g2hjZZ,g2hjZga,g2hjgaga,g2hjgg,g2hjggZ,&
&    g2hjhiZ, BR_hjinvisible,BR_hjhihi)

!-Run HS on the original experimental data in order to evaluate the model predictions	  
  call run_HiggsSignals(1, Chisq_mu, Chisq_mh, Chisq, ndf, Pvalue)

!-Get the number of the peak-observables (Don't care about ntotal, npeakmh, nmpred, nanalyses)
  call get_number_of_observables(ntotal, npeakmu, npeakmh, nmpred, nanalyses)

! We now want to set different efficiencies as evaluated e.g. from a MC analysis.
! Note that these are given as relative changes with respect to the SM efficiencies, i.e.
! in the MC you should both test your model and the SM with the analysis cuts and
! estimate the ratio eta = efficiency(model)/efficiency(SM) for each Higgs channel.
!
! As an example, we want to give all VBF,WH,ZH channels a factor of etaZ in all implemented
! searches where these channels are present.

 do ll=1,101
 etaZ = 0.0D0 + (ll-1)*0.025D0
!-Loop over the number of peak observables
  do kk=1,npeakmu
!--Get for each peak observable its unique ID:
   call get_ID_of_peakobservable(kk, ID)
!--Get the predicted signal strength modifier (mupred) for this peak observable:
!   call get_peakinfo_from_HSresults(ID, mupred, domH, nHcomb)
!--Get channel information
   call get_peak_channels(ID, Nc, channel_IDs, efficiencies)
!   write(*,*) ID, Nc
!   write(*,*) channel_IDs
!   write(*,*) efficiencies
   allocate(eff_ratio(Nc))
   eff_ratio = 1.0D0
   do jj=lbound(channel_IDs,dim=1),ubound(channel_IDs,dim=1)
!---Get the ID for the production mode (2:VBF, 3:WH, 4:ZH):   
    pID = int((channel_IDs(jj)-modulo(channel_IDs(jj),10))/dble(10))
    if(pID.eq.2.or.pID.eq.3.or.pID.eq.4) then
     eff_ratio(jj)=etaZ
    endif
   enddo 
!---Hand over the efficiency ratios:
   call assign_modelefficiencies_to_peak(ID, Nc, eff_ratio) 
   deallocate(channel_IDs,efficiencies,eff_ratio)
  enddo
    
!-Set the HiggsSignals input (again!)
  call HiggsBounds_neutral_input_effC(Mh,GammaTotal,&
&    g2hjss_s,g2hjss_p,g2hjcc_s,g2hjcc_p,g2hjbb_s,g2hjbb_p,&
&    g2hjtt_s,g2hjtt_p,&
&    g2hjmumu_s,g2hjmumu_p,g2hjtautau_s,g2hjtautau_p,&
&    g2hjWW,g2hjZZ,g2hjZga,g2hjgaga,g2hjgg,g2hjggZ,&
&    g2hjhiZ, BR_hjinvisible,BR_hjhihi)

  call run_HiggsSignals(1, Chisq_mu, Chisq_mh, Chisq, ndf, Pvalue)

   mu_average = 0.0D0
  do kk=1,npeakmu
!--Get for each peak observable its unique ID:
   call get_ID_of_peakobservable(kk, ID)
!--Get the predicted signal strength modifier (mupred) for this peak observable:
   call get_peakinfo_from_HSresults(ID, mupred, domH, nHcomb)
   mu_average = mu_average + mupred
  enddo
  mu_average = mu_average / npeakmu

  write(21,*) g2hjWW, etaZ, Chisq_mu, mu_average

 enddo
 close(21)
 enddo


 call finish_HiggsSignals
 
!--------------------------------------------------------------------------------------
end program HS_efficiencies
!--------------------------------------------------------------------------------------
