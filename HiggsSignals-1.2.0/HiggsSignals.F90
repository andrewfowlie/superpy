!--------------------------------------------------------------------
! This file is part of HiggsSignals (TS 31/01/2013)
!--------------------------------------------------------------------
program HiggsSignals
!--------------------------------------------------------------------
! Creates the command-line executable of HiggsSignals 
!
! HiggsSignals performs a chi^2 test of Higgs signal rate and mass 
! predictions of an arbitrary Higgs sector with measurements from 
! the Tevatron and LHC experiments.
!--------------------------------------------------------------------
 use usefulbits, only : inputmethod,np,Hneut,Hplus,whichanalyses
 use input, only : do_input
 use usefulbits_hs, only : output_level, HiggsSignals_info, runmode, Exptdir
 use io, only : setup_input_for_hs, do_output_for_hs

 implicit none

 double precision :: Pvalue, Chisq, Chisq_mu, Chisq_mh
 integer :: ndf, mode

 inputmethod='datfile'
 output_level=0
 whichanalyses='onlyH'
  
 call setup_input_for_hs
 
 call initialize_HiggsSignals(np(Hneut),np(Hplus),Exptdir)
 
 call do_input
 
 if(runmode.eq.'peak') then
  mode = 1
 elseif(runmode.eq.'mass') then
  mode = 2
 elseif(runmode.eq.'both') then
  mode = 3
 else
  stop"Error: runmode is unknown. Please specify as peak, mass or both." 
 endif
 
 call run_HiggsSignals( mode, Chisq_mu, Chisq_mh, Chisq, ndf, Pvalue)
 
 call do_output_for_hs
 
 call finish_HiggsSignals

end program