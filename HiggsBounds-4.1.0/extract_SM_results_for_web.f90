program extract_SM_results_for_web
! compile with the command:
! gfortran extract_SM_results_for_web.f90 -o extract_SM_results_for_web -L. -lHB

 implicit none
 ! SM functions from the HiggsBounds library
 double precision :: SMBR_HWW
 double precision :: SMBR_HZZ
 double precision :: SMBR_Hbb
 double precision :: SMBR_Htautau
 double precision :: SMBR_Hgamgam
 double precision :: SMBR_Hgg
 double precision :: SMGamma_h
 double precision :: SMGamma_tWpb

 double precision :: SMBR_HZgam
 double precision :: SMBR_Htoptop
 double precision :: SMBR_Hcc
 double precision :: SMBR_Hss
 double precision :: SMBR_Hmumu

 double precision :: SMCS_tev_HW
 double precision :: SMCS_tev_HZ
 double precision :: SMCS_tev_gg_H
 double precision :: SMCS_tev_bb_H
 double precision :: SMCS_tev_vbf_H
 double precision :: SMCS_tev_bg_Hb
 double precision :: SMCS_tev_ttH

 double precision :: SMCS_lhc7_HW
 double precision :: SMCS_lhc7_HZ
 double precision :: SMCS_lhc7_gg_H
 !double precision :: SMCS_lhc7_bb_H !not in code yet
 double precision :: SMCS_lhc7_vbf_H
 !double precision :: SMCS_lhc7_bg_Hb !not in code yet
 double precision :: SMCS_lhc7_ttH

 double precision :: SMCS_lhc8_HW
 double precision :: SMCS_lhc8_HZ
 double precision :: SMCS_lhc8_gg_H
 !double precision :: SMCS_lhc8_bb_H !not in code yet
 double precision :: SMCS_lhc8_vbf_H
 !double precision :: SMCS_lhc8_bg_Hb !not in code yet
 double precision :: SMCS_lhc8_ttH


 ! internal
 double precision :: Mhi,mt
 integer :: number_args 
 character(LEN=100) :: temp
 character(LEN=9) :: f_out
 double precision :: SMtemp
 integer :: iargc

 number_args = IARGC() 
   
 if(number_args.ne.1)then
  stop 'program extract_SM_results_for_web has wrong number of arguments'
 endif

 temp=""
 call GETARG(1,temp)
 read(temp,*) Mhi

 mt=173.1D0

 ! call subroutine initialize_HiggsBounds with dummy arguments
 call initialize_HiggsBounds(1,0,'LandH')

 write(*,*)''
 write(*,*)'~ some of the Standard Model quantities used internally in HiggsBounds ~'
 write(*,*)''
 write(*,*)'For a SM Higgs with a mass (in GeV) of ',Mhi
 write(*,*)'the SM branching ratios are'

 f_out='(a,G16.6)'

 SMtemp=              SMBR_HWW(Mhi)     ; call flush(6)
 write(*,fmt=f_out)'  SMBR_HWW      = ' , SMtemp
 SMtemp=              SMBR_HZZ(Mhi)     ; call flush(6)
 write(*,fmt=f_out)'  SMBR_HZZ      = ' , SMtemp
 SMtemp=              SMBR_Hbb(Mhi)     ; call flush(6)
 write(*,fmt=f_out)'  SMBR_Hbb      = ' , SMtemp
 SMtemp=              SMBR_Htautau(Mhi) ; call flush(6)
 write(*,fmt=f_out)'  SMBR_Htautau  = ' , SMtemp
 SMtemp=              SMBR_Hgamgam(Mhi) ; call flush(6)
 write(*,fmt=f_out)'  SMBR_Hgamgam  = ' , SMtemp
 SMtemp=              SMBR_Hgg(Mhi)     ; call flush(6)
 write(*,fmt=f_out)'  SMBR_Hgg      = ' , SMtemp

 SMtemp=              SMBR_HZgam(Mhi)     ; call flush(6)
 write(*,fmt=f_out)'  SMBR_HZgam    = ' , SMtemp
 SMtemp=              SMBR_HZZ(Mhi)     ; call flush(6)
 write(*,fmt=f_out)'  SMBR_Htoptop  = ' , SMtemp
 SMtemp=              SMBR_Hcc(Mhi)     ; call flush(6)
 write(*,fmt=f_out)'  SMBR_Hcc      = ' , SMtemp
 SMtemp=              SMBR_Hss(Mhi)     ; call flush(6)
 write(*,fmt=f_out)'  SMBR_Hss      = ' , SMtemp
 SMtemp=              SMBR_Hmumu(Mhi)     ; call flush(6)
 write(*,fmt=f_out)'  SMBR_Hmumu    = ' , SMtemp

 write(*,*)'the total SM decay width (in Gev) is'
 SMtemp=              SMGamma_h(Mhi)    ; call flush(6)
 write(*,fmt=f_out)'  SMGamma_h     = ' , SMtemp ; call flush(6)

 write(*,*)'the SM top decay width (in Gev) into a W boson and b quark is'
 SMtemp=              SMGamma_tWpb(mt)    ; call flush(6)
 write(*,fmt=f_out)'  SMGamma_tWpb     = ' , SMtemp ; call flush(6)
 write(*,*)'     for a top mass of ',mt

 write(*,*)'the SM Tevatron hadronic cross sections (in pb) are'

 SMtemp=              SMCS_tev_HW(Mhi)    ; call flush(6)
 write(*,fmt=f_out)'  SMCS_tev_HW        = ' , SMtemp
 SMtemp=              SMCS_tev_HZ(Mhi)    ; call flush(6)
 write(*,fmt=f_out)'  SMCS_tev_HZ        = ' , SMtemp
 SMtemp=              SMCS_tev_gg_H(Mhi)    ; call flush(6)
 write(*,fmt=f_out)'  SMCS_tev_gg_H      = ' , SMtemp
 SMtemp=              SMCS_tev_bb_H(Mhi)    ; call flush(6)
 write(*,fmt=f_out)'  SMCS_tev_bb_H      = ' , SMtemp
 SMtemp=              SMCS_tev_vbf_H(Mhi)    ; call flush(6)
 write(*,fmt=f_out)'  SMCS_tev_vbf_H     = ' , SMtemp
 SMtemp=              SMCS_tev_bg_Hb(Mhi)    ; call flush(6)
 write(*,fmt=f_out)'  SMCS_tev_bg_Hb     = ' , SMtemp
 SMtemp=              SMCS_tev_ttH(Mhi)      ; call flush(6)
 write(*,fmt=f_out)'  SMCS_tev_ttH       = ' , SMtemp

 write(*,*)'the SM LHC hadronic cross sections (in pb) at 7TeV are'

 SMtemp=              SMCS_lhc7_HW(Mhi)    ; call flush(6)
 write(*,fmt=f_out)'  SMCS_lhc7_HW       = ' , SMtemp
 SMtemp=              SMCS_lhc7_HZ(Mhi)    ; call flush(6)
 write(*,fmt=f_out)'  SMCS_lhc7_HZ       = ' , SMtemp
 SMtemp=              SMCS_lhc7_gg_H(Mhi)    ; call flush(6)
 write(*,fmt=f_out)'  SMCS_lhc7_gg_H     = ' , SMtemp
 !SMtemp=              SMCS_lhc7_bb_H(Mhi)    ; call flush(6)
 !write(*,fmt=f_out)'  SMCS_lhc7_bb_H     = ' , SMtemp
 SMtemp=              SMCS_lhc7_vbf_H(Mhi)    ; call flush(6)
 write(*,fmt=f_out)'  SMCS_lhc7_vbf_H    = ' , SMtemp
 !SMtemp=              SMCS_lhc7_bg_Hb(Mhi)    ; call flush(6)
 !write(*,fmt=f_out)'  SMCS_lhc7_bg_Hb    = ' , SMtemp
 SMtemp=              SMCS_lhc7_ttH(Mhi)      ; call flush(6)
 write(*,fmt=f_out)'  SMCS_lhc7_ttH      = ' , SMtemp


 write(*,*)'the SM LHC hadronic cross sections (in pb) at 8TeV are'

 SMtemp=              SMCS_lhc8_HW(Mhi)    ; call flush(6)
 write(*,fmt=f_out)'  SMCS_lhc8_HW       = ' , SMtemp
 SMtemp=              SMCS_lhc8_HZ(Mhi)    ; call flush(6)
 write(*,fmt=f_out)'  SMCS_lhc8_HZ       = ' , SMtemp
 SMtemp=              SMCS_lhc8_gg_H(Mhi)    ; call flush(6)
 write(*,fmt=f_out)'  SMCS_lhc8_gg_H     = ' , SMtemp
 !SMtemp=              SMCS_lhc7_bb_H(Mhi)    ; call flush(6)
 !write(*,fmt=f_out)'  SMCS_lhc7_bb_H     = ' , SMtemp
 SMtemp=              SMCS_lhc8_vbf_H(Mhi)    ; call flush(6)
 write(*,fmt=f_out)'  SMCS_lhc8_vbf_H    = ' , SMtemp
 !SMtemp=              SMCS_lhc7_bg_Hb(Mhi)    ; call flush(6)
 !write(*,fmt=f_out)'  SMCS_lhc7_bg_Hb    = ' , SMtemp
 SMtemp=              SMCS_lhc8_ttH(Mhi)      ; call flush(6)
 write(*,fmt=f_out)'  SMCS_lhc8_ttH      = ' , SMtemp


 !write(*,*)
 !write(*,*)'These results were obtained from' 
 !write(*,*)'the program HDECAY (arXiv:hep-ph/9704448)'
 !write(*,*)'and cross sections collated by the'
 !write(*,*)'and the TeV4LHC Higgs Working Group(arXiv:hep-ph/0612172)'
 write(*,*)
 write(*,*)'See HiggsBounds manual for definitions and references.'
 write(*,*)
 write(*,*)'Note that, for the majority of SM quantities used internally by HiggsBounds,'
 write(*,*)'we make use of existing calculations, including the program HDECAY (arXiv:hep-ph/9704448),'
 write(*,*)'those collated by the TeV4LHC Higgs Working Group (arXiv:hep-ph/0612172)'
 write(*,*)'and those collated by the LHC Higgs Cross Section Working Group (arXiv:1101.0593).'
 write(*,*)
 write(*,*)'For a full list of the appropriate references in .bib and .bbl format, see'
 write(*,*)'http://projects.hepforge.org/higgsbounds/downloads.html'


 call finish_HiggsBounds

end program extract_SM_results_for_web
