! Contains the functions which allow the user to access the Standard Model Higgs
! Branching ratios, total decay width and cross sections.
!
! This file is part of HiggsBounds
!  -KW
!*********************************
function SMBR_HWW(Mh)
  use theory_BRfunctions
  implicit none
  double precision,intent(in) :: Mh
  double precision :: SMBR_HWW

  call testBRSM(Mh)
  
  SMBR_HWW=BRSM_HWW(Mh,.True.)

end function SMBR_HWW
!*********************************
function SMBR_HZZ(Mh)
  use theory_BRfunctions
  implicit none
  double precision,intent(in) :: Mh
  double precision :: SMBR_HZZ

  call testBRSM(Mh)
  
  SMBR_HZZ=BRSM_HZZ(Mh,.True.)

end function SMBR_HZZ
!*********************************
function SMBR_Hbb(Mh)
  use theory_BRfunctions
  implicit none
  double precision,intent(in) :: Mh
  double precision :: SMBR_Hbb
 
  call testBRSM(Mh)
  
  SMBR_Hbb=BRSM_Hbb(Mh,.True.)

end function SMBR_Hbb
!*********************************
function SMBR_Htautau(Mh)
  use theory_BRfunctions
  implicit none
  double precision,intent(in) :: Mh
  double precision :: SMBR_Htautau

  call testBRSM(Mh)
  
  SMBR_Htautau=BRSM_Htautau(Mh,.True.)
 
end function SMBR_Htautau
!*********************************
function SMBR_Hgamgam(Mh)
  use theory_BRfunctions
  implicit none
  double precision,intent(in) :: Mh
  double precision :: SMBR_Hgamgam

  call testBRSM(Mh)
  
  SMBR_Hgamgam=BRSM_Hgaga(Mh,.True.)

end function SMBR_Hgamgam
!*********************************
function SMBR_Hgg(Mh)
  use theory_BRfunctions
  implicit none
  double precision,intent(in) :: Mh
  double precision :: SMBR_Hgg

  call testBRSM(Mh)
  
  SMBR_Hgg=BRSM_Hgg(Mh,.True.)
end function SMBR_Hgg
!*********************************
function SMBR_Htoptop(Mh)
  use theory_BRfunctions
  implicit none
  double precision,intent(in) :: Mh
  double precision :: SMBR_Htoptop

  call testBRSM(Mh)
  
  SMBR_Htoptop=BRSM_Htoptop(Mh,.True.) 
end function SMBR_Htoptop
!*********************************
function SMBR_Hcc(Mh)
  use theory_BRfunctions
  implicit none
  double precision,intent(in) :: Mh
  double precision :: SMBR_Hcc

  call testBRSM(Mh)
  
  SMBR_Hcc=BRSM_Hcc(Mh,.True.)
 end function SMBR_Hcc
!*********************************
function SMBR_Hss(Mh)
  use theory_BRfunctions
  implicit none
  double precision,intent(in) :: Mh
  double precision :: SMBR_Hss

  call testBRSM(Mh)
  
  SMBR_Hss=BRSM_Hss(Mh,.True.)
end function SMBR_Hss
!*********************************
function SMBR_Hmumu(Mh)
  use theory_BRfunctions
  implicit none
  double precision,intent(in) :: Mh
  double precision :: SMBR_Hmumu

  call testBRSM(Mh)
  
  SMBR_Hmumu=BRSM_Hmumu(Mh,.True.)
end function SMBR_Hmumu
!*********************************
function SMBR_HZgam(Mh)
  use theory_BRfunctions
  implicit none
  double precision,intent(in) :: Mh
  double precision :: SMBR_HZgam
  
  call testBRSM(Mh)
  
  SMBR_HZgam=BRSM_HZga(Mh,.True.)
end function SMBR_HZgam
!*********************************
function SMGamma_h(Mh)
  use theory_BRfunctions
  implicit none
  double precision,intent(in) :: Mh
  double precision :: SMGamma_h

  call testBRSM(Mh)
  
  SMGamma_h=BRSM_GammaTot(Mh,.True.)

end function SMGamma_h
!*********************************
function SMGamma_tWpb(mt)
  use theory_BRfunctions
  implicit none
  double precision,intent(in) :: mt
  double precision :: SMGamma_tWpb

  call testBRSM(mt)

  SMGamma_tWpb=BRSM_Gamma_tWpb(mt,.True.)

end function SMGamma_tWpb
!*********************************
function SMCS_tev_HW(Mh)
!note: in pb
 use theory_XS_SM_functions
  implicit none
  double precision,intent(in) :: Mh
  double precision :: SMCS_tev_HW
  
  call testBRSM(Mh)
  
  SMCS_tev_HW=1.0D-3*XS_tev_HW_SM(Mh,.True.)
end function SMCS_tev_HW
!*********************************
function SMCS_tev_HZ(Mh)
!note: in pb
 use theory_XS_SM_functions
  implicit none
  double precision,intent(in) :: Mh
  double precision :: SMCS_tev_HZ
  
  call testBRSM(Mh)
  
  SMCS_tev_HZ=1.0D-3*XS_tev_HZ_SM(Mh,.True.)
end function SMCS_tev_HZ
!*********************************
function SMCS_tev_gg_H(Mh)
!note: in pb
 use theory_XS_SM_functions
  implicit none
  double precision,intent(in) :: Mh
  double precision :: SMCS_tev_gg_H
  
  call testBRSM(Mh)
  
!  SMCS_tev_gg_H=1.0D-3*XS_tev_gg_H_SM(Mh) !use this for outdated CS
  SMCS_tev_gg_H=1.0D-3*XS_tev_gg_H_SM(Mh,.True.) !use this for updated (31/03/2011) CS
end function SMCS_tev_gg_H
!*********************************
function SMCS_tev_bb_H(Mh)
!note: in pb
 use theory_XS_SM_functions
  implicit none
  double precision,intent(in) :: Mh
  double precision :: SMCS_tev_bb_H
  
  call testBRSM(Mh)
  
  SMCS_tev_bb_H=1.0D-3*XS_tev_bb_H_SM(Mh,.True.)
end function SMCS_tev_bb_H
!*********************************
function SMCS_tev_vbf_H(Mh)
!note: in pb
 use theory_XS_SM_functions
  implicit none
  double precision,intent(in) :: Mh
  double precision :: SMCS_tev_vbf_H
  
  call testBRSM(Mh)
  
  SMCS_tev_vbf_H=1.0D-3*XS_tev_vbf_SM(Mh,.True.)
end function SMCS_tev_vbf_H
!*********************************
function SMCS_tev_bg_Hb(Mh) 
!note: in pb
 use theory_XS_SM_functions
  implicit none
  double precision,intent(in) :: Mh
  double precision :: SMCS_tev_bg_Hb
  
  call testBRSM(Mh)
  
  SMCS_tev_bg_Hb=1.0D-3*XS_tev_bg_Hb_SM(Mh,.True.)
end function SMCS_tev_bg_Hb
!*********************************
function SMCS_tev_ttH(Mh) 
!note: in pb
 use theory_XS_SM_functions
  implicit none
  double precision,intent(in) :: Mh
  double precision :: SMCS_tev_ttH
  
  call testBRSM(Mh)
  
  SMCS_tev_ttH=1.0D-3*XS_tev_ttH_SM(Mh,.True.)
end function SMCS_tev_ttH
!*********************************
function SMCS_lhc7_HW(Mh) 
!note: in pb
 use theory_XS_SM_functions
  implicit none
  double precision,intent(in) :: Mh
  double precision :: SMCS_lhc7_HW
  
  call testBRSM(Mh)
  
  SMCS_lhc7_HW=XS_lhc7_HW_SM(Mh,.True.)
end function SMCS_lhc7_HW
!*********************************
function SMCS_lhc7_HZ(Mh) 
!note: in pb
 use theory_XS_SM_functions
  implicit none
  double precision,intent(in) :: Mh
  double precision :: SMCS_lhc7_HZ
  
  call testBRSM(Mh)
  
  SMCS_lhc7_HZ=XS_lhc7_HZ_SM(Mh,.True.)
end function SMCS_lhc7_HZ
!*********************************
function SMCS_lhc7_gg_H(Mh) 
!note: in pb
 use theory_XS_SM_functions
  implicit none
  double precision,intent(in) :: Mh
  double precision :: SMCS_lhc7_gg_H
  
  call testBRSM(Mh)
  
  SMCS_lhc7_gg_H=XS_lhc7_gg_H_SM(Mh,.True.)
end function SMCS_lhc7_gg_H
!*********************************
function SMCS_lhc7_bb_H(Mh) 
!note: in pb
 use theory_XS_SM_functions
  implicit none
  double precision,intent(in) :: Mh
  double precision :: SMCS_lhc7_bb_H
  
  call testBRSM(Mh)
  
  SMCS_lhc7_bb_H=XS_lhc7_bb_H_SM(Mh,.True.)
end function SMCS_lhc7_bb_H
!*********************************
function SMCS_lhc7_vbf_H(Mh) 
!note: in pb
 use theory_XS_SM_functions
  implicit none
  double precision,intent(in) :: Mh
  double precision :: SMCS_lhc7_vbf_H
  
  call testBRSM(Mh)
  
  SMCS_lhc7_vbf_H=XS_lhc7_vbf_SM(Mh,.True.)
end function SMCS_lhc7_vbf_H
!*********************************
function SMCS_lhc7_ttH(Mh) 
!note: in pb
 use theory_XS_SM_functions
  implicit none
  double precision,intent(in) :: Mh
  double precision :: SMCS_lhc7_ttH
  
  call testBRSM(Mh)
  
  SMCS_lhc7_ttH=XS_lhc7_ttH_SM(Mh,.True.)
end function SMCS_lhc7_ttH
!*********************************

function SMCS_lhc8_HW(Mh) 
!note: in pb
 use theory_XS_SM_functions
  implicit none
  double precision,intent(in) :: Mh
  double precision :: SMCS_lhc8_HW
  
  call testBRSM(Mh)
  
  SMCS_lhc8_HW=XS_lhc8_HW_SM(Mh,.True.)
end function SMCS_lhc8_HW
!*********************************
function SMCS_lhc8_HZ(Mh) 
!note: in pb
 use theory_XS_SM_functions
  implicit none
  double precision,intent(in) :: Mh
  double precision :: SMCS_lhc8_HZ
  
  call testBRSM(Mh)
  
  SMCS_lhc8_HZ=XS_lhc8_HZ_SM(Mh,.True.)
end function SMCS_lhc8_HZ
!*********************************
function SMCS_lhc8_gg_H(Mh) 
!note: in pb
 use theory_XS_SM_functions
  implicit none
  double precision,intent(in) :: Mh
  double precision :: SMCS_lhc8_gg_H
  
  call testBRSM(Mh)
  
  SMCS_lhc8_gg_H=XS_lhc8_gg_H_SM(Mh,.True.)
end function SMCS_lhc8_gg_H
!*********************************
function SMCS_lhc8_bb_H(Mh) 
!note: in pb
 use theory_XS_SM_functions
  implicit none
  double precision,intent(in) :: Mh
  double precision :: SMCS_lhc8_bb_H
  
  call testBRSM(Mh)
  
  SMCS_lhc8_bb_H=XS_lhc8_bb_H_SM(Mh,.True.)
end function SMCS_lhc8_bb_H
!*********************************
function SMCS_lhc8_vbf_H(Mh) 
!note: in pb
 use theory_XS_SM_functions
  implicit none
  double precision,intent(in) :: Mh
  double precision :: SMCS_lhc8_vbf_H
  
  call testBRSM(Mh)
  
  SMCS_lhc8_vbf_H=XS_lhc8_vbf_SM(Mh,.True.)
end function SMCS_lhc8_vbf_H
!*********************************
function SMCS_lhc8_ttH(Mh) 
!note: in pb
 use theory_XS_SM_functions
  implicit none
  double precision,intent(in) :: Mh
  double precision :: SMCS_lhc8_ttH
  
  call testBRSM(Mh)
  
  SMCS_lhc8_ttH=XS_lhc8_ttH_SM(Mh,.True.)
end function SMCS_lhc8_ttH
!*********************************

subroutine testBRSM(M)
  use theory_BRfunctions
  implicit none
  double precision,intent(in) :: M

  if(.not.allocated(BRSM))then
   write(*,*)'You can only use this function between'
   write(*,*)'calling the subroutines'
   write(*,*)'initialize_HiggsBounds and'
   write(*,*)'finish_HiggsBounds'
   stop'error (see standard output for more info)'
  endif

  !if(M.gt.0.0D0)then
   !ok
  !else !M is negative or NaN
  ! stop'wrong mass given to function'
  !endif
end subroutine testBRSM
!*********************************
! OLD FUNCTIONS
!*********************************
function SMCS_tev_qq_HZ(Mh) 
  implicit none
  double precision,intent(in) :: Mh
  double precision :: SMCS_tev_qq_HZ
  double precision :: SMCS_tev_HZ

  write(*,*)'Note: function SMCS_tev_qq_HZ has been'
  write(*,*)'superceded by function SMCS_tev_HZ'
  SMCS_tev_qq_HZ=SMCS_tev_HZ(Mh)

end function SMCS_tev_qq_HZ
!*********************************
function SMCS_tev_qq_HW(Mh) 
  implicit none
  double precision,intent(in) :: Mh
  double precision :: SMCS_tev_qq_HW
  double precision :: SMCS_tev_HW

  write(*,*)'Note: function SMCS_tev_qq_HW has been'
  write(*,*)'superceded by function SMCS_tev_HW'
  SMCS_tev_qq_HW=SMCS_tev_HW(Mh)

end function SMCS_tev_qq_HW
!*********************************
