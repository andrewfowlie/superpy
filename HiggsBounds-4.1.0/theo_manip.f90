! This file is part of HiggsBounds
!  -KW
!******************************************************************
module theo_manip
!******************************************************************
 !use S95tables_type1
 use usefulbits, only : ndat,np,Hneut,Hplus,theo,partR,hadroncolliderextras,pdesc
 
 implicit none

 type(hadroncolliderextras) :: tevS(1)
 type(hadroncolliderextras) :: lhc7S(1)
 type(hadroncolliderextras) :: lhc8S(1)

 contains

 !******************************************************************
 subroutine complete_theo
 !decides what has to be done to the input and calls the appropriate 
 !subroutines
 !******************************************************************
  use usefulbits, only : whichanalyses,whichinput,ndat
  implicit none
  !--------------------------------------internal
  integer :: i,j,jj
  !----------------------------------------------

  if(np(Hneut)>0)then !none if this is needed for the charged Higgs sector yet

   select case(whichinput)        
   case('effC')  
    call csratios_from_g2
    call cp_from_g2
    call br_from_g2
   case('SLHA')
    call csratios_from_g2
    call cp_from_g2
   case('hadr','part')
   case default
    stop'error in subroutine complete_theo (1)' 
   end select               
    
  
   do jj=1,ndat     ! filling the other half of XS_hjhi_ratio
    do j=2,np(Hneut)
     do i=1,j-1          
      theo(jj)%lep%XS_hjhi_ratio(i,j) = theo(jj)%lep%XS_hjhi_ratio(j,i)               
     enddo
    enddo
   enddo  
  endif

  call check_dataset !involves the charged Higgs sector
 
  if(np(Hneut)>0)then !none if this is needed for the charged Higgs sector yet
   select case(whichanalyses)
   case('onlyH','LandH','onlyP','list ') ! everything which involves Tevatron and LHC tables

    call fill_theo_SM ! n.b. there's no LEP SM cross sections at the moment
  
     select case(whichinput)
     case('part','effC','SLHA') ! everything except option 'hadr', where had XS ratios are inputted directly
      call XS_from_partR
     case('hadr')
     case default
      stop'error in subroutine complete_theo (2)'  
     end select
   case('onlyL')
   case default       
    stop'error in subroutine complete_theo (3)'             
   end select   
  endif

 end subroutine complete_theo
 !******************************************************************
 subroutine recalculate_theo_for_datapoint(n)
 ! Does the same as complete_theo but just for the datapoint n.
  use usefulbits, only : whichanalyses,whichinput
  implicit none
  integer, intent(in) :: n

  if(np(Hneut)>0) then !none if this is needed for the charged Higgs sector yet
   select case(whichinput)
   case('effC')
    call csratios_from_g2_for_datapoint(n)
    call br_from_g2_for_datapoint(n)
   case('SLHA')
    call csratios_from_g2_for_datapoint(n)
   case('hadr','part')
   case default
    stop'error in subroutine recalculate_theo_for_datapoint (1)' 
   end select               
  endif
    
  call check_dataset !involves the charged Higgs sector
 
  if(np(Hneut)>0)then !none if this is needed for the charged Higgs sector yet
   select case(whichanalyses)
   case('onlyH','LandH','onlyP','list ') ! everything which involves Tevatron and LHC tables
    call fill_theo_SM_for_datapoint(n) ! n.b. there's no LEP SM cross sections at the moment
    select case(whichinput)
     case('part','effC','SLHA') ! everything except option 'hadr', where had XS ratios are inputted directly
      call XS_from_partR_for_datapoint(n)
     case('hadr')
     case default
      stop'error in subroutine recalculate_theo_for_datapoint (2)'  
     end select
   case('onlyL')
   case default       
    stop'error in subroutine recalculate_theo_for_datapoint (3)'             
   end select   
  endif

 end subroutine recalculate_theo_for_datapoint
 !****************************************************************** 
 subroutine csratios_from_g2 
 ! uses the effective couplings contained in g2 to calculate 
 ! partonic cross section ratios, some hadronic cross section ratios
 !*****************************************************************
  use usefulbits, only : ndat
  implicit none
  !--------------------------------------internal
  integer :: jj
  !---------------------------------------------

  if(np(Hneut)<1)stop'error in csratios_from_g2  (np(Hneut))'

   do jj=1,ndat
    call csratios_from_g2_for_datapoint(jj)
   enddo

 end subroutine csratios_from_g2
 !****************************************************************** 
 subroutine csratios_from_g2_for_datapoint(jj) 
 ! uses the effective couplings contained in g2 to calculate 
 ! partonic cross section ratios, some hadronic cross section ratios
 !*****************************************************************
  use usefulbits, only : g2
  use theory_colliderSfunctions
  use S95tables, only : inrange
  implicit none
  integer, intent(in) :: jj
  !--------------------------------------internal
  integer :: i
  double precision :: TEVSM_ZZ_contrib_to_VBF,TEVSM_WW_contrib_to_VBF
  double precision :: Mhi
  !---------------------------------------------

! relative contributuion of WW- and ZZ-fusion to VBF (in LO) for
! p p-bar collisions at SqrtS=1.96 TeV (calcuated by T. Figy with VBFNLO):s
   TEVSM_ZZ_contrib_to_VBF=0.23D0
   TEVSM_WW_contrib_to_VBF=0.77D0

   do i=1,np(Hneut)
                
    theo(jj)%lep%XS_hjZ_ratio(i)      =  g2(jj)%hjZZ(i) 
    theo(jj)%lep%XS_bbhj_ratio(i)     =  g2(jj)%hjbb_s(i)+g2(jj)%hjbb_p(i)!nb tables at the moment can not be applied to mixed CP Higgs
    theo(jj)%lep%XS_tautauhj_ratio(i) =  g2(jj)%hjtautau_s(i)+g2(jj)%hjtautau_p(i)!nb tables at the moment can not be applied to mixed CP Higgs

    partR(jj)%bg_hjb(i)             =   g2(jj)%hjbb_s(i)+g2(jj)%hjbb_p(i)
    theo(jj)%tev%XS_vbf_ratio(i)   =   g2(jj)%hjWW(i)*TEVSM_WW_contrib_to_VBF &        
                               &       + g2(jj)%hjZZ(i)*TEVSM_ZZ_contrib_to_VBF
    theo(jj)%tev%XS_tthj_ratio(i)  =   g2(jj)%hjtoptop_s(i)+g2(jj)%hjtoptop_p(i) !nb tev tables at the moment can only use CP even Higgs

    Mhi=theo(jj)%particle(Hneut)%M(i)

    if(inrange(Mhi,'LHC7'))then
     theo(jj)%lhc7%XS_vbf_ratio(i)   =   g2(jj)%hjWW(i)*lhc7_rHVBF_WW(Mhi) &        
                                 &        + g2(jj)%hjZZ(i)*lhc7_rHVBF_ZZ(Mhi)
    else 
     theo(jj)%lhc7%XS_vbf_ratio(i)   = 0.0D0  
    endif

    theo(jj)%lhc7%XS_tthj_ratio(i)  =   g2(jj)%hjtoptop_s(i)+g2(jj)%hjtoptop_p(i) !nb no tables need this at the moment

!	  We are using 7 TeV ratios for VBF contribution from WW/ZZ at the moment also 
!     for LHC 8 TeV cross sections
    if(inrange(Mhi,'LHC8'))then
     theo(jj)%lhc8%XS_vbf_ratio(i)   =   g2(jj)%hjWW(i)*lhc7_rHVBF_WW(Mhi) &        
                              &        + g2(jj)%hjZZ(i)*lhc7_rHVBF_ZZ(Mhi)
    else 
     theo(jj)%lhc8%XS_vbf_ratio(i)   = 0.0D0  
    endif

    theo(jj)%lhc8%XS_tthj_ratio(i)  =   g2(jj)%hjtoptop_s(i)+g2(jj)%hjtoptop_p(i) !nb no tables need this at the moment
           
    partR(jj)%qq_hjWp(:,i)      = g2(jj)%hjWW(i)
    partR(jj)%qq_hjWm(:,i)      = g2(jj)%hjWW(i) 
    partR(jj)%gg_hj(i)          = g2(jj)%hjgg(i)
    partR(jj)%qq_hj(5,i)        = g2(jj)%hjbb_s(i)+g2(jj)%hjbb_p(i)        
    partR(jj)%qq_hjZ(:,i)       = g2(jj)%hjZZ(i)   
    partR(jj)%gg_hjZ(i)         = g2(jj)%hjggZ(i)        
        
   enddo 

   theo(jj)%lep%XS_hjhi_ratio=g2(jj)%hjhiZ! note only half of XS_hjhi_ratio is filled here

 end subroutine csratios_from_g2_for_datapoint
 !****************************************************************** 
 subroutine cp_from_g2 
 ! uses the effective couplings contained in g2 to calculate 
 ! cp of neutral higgs
 !*****************************************************************
  use usefulbits, only : g2,ndat,vsmall
  implicit none
  !--------------------------------------internal
  integer :: i,jj
  double precision :: max_hjff_s,max_hjff_p
  !---------------------------------------------

  if(np(Hneut)<1)stop'error in cp_from_g2  (np(Hneut))'

   do jj=1,ndat
 
    do i=1,np(Hneut)
     max_hjff_s=max(g2(jj)%hjss_s(i),g2(jj)%hjcc_s(i),g2(jj)%hjbb_s(i), &
         &   g2(jj)%hjtoptop_s(i),g2(jj)%hjmumu_s(i),g2(jj)%hjtautau_s(i))

     max_hjff_p=max(g2(jj)%hjss_p(i),g2(jj)%hjcc_p(i),g2(jj)%hjbb_p(i), &
         &   g2(jj)%hjtoptop_p(i),g2(jj)%hjmumu_p(i),g2(jj)%hjtautau_p(i))

     if(     max_hjff_p .lt. vsmall )then !CP even
      theo(jj)%CP_value(i) =  1
     elseif( max_hjff_s .lt. vsmall )then !CP odd
      theo(jj)%CP_value(i) = -1
     else                              !mixed CP
      theo(jj)%CP_value(i) =  0
     endif
    enddo 
 
   enddo

 end subroutine cp_from_g2

 !****************************************************************** 
 subroutine br_from_g2 
 ! uses the effective couplings contained in g2 to calculate 
 ! branching ratios
 !*****************************************************************
  use usefulbits, only : np,Hneut,ndat
  implicit none
  !--------------------------------------internal
  integer :: jj
  !---------------------------------------------

  if(np(Hneut)<1)stop'error in br_from_g2 (np(Hneut))'

   do jj=1,ndat
    call br_from_g2_for_datapoint(jj)
   enddo

 end subroutine br_from_g2
 !*****************************************************************
 subroutine br_from_g2_for_datapoint(jj)
 ! uses the effective couplings contained in g2 to calculate 
 ! branching ratios
 !*****************************************************************
  use theory_BRfunctions
  use S95tables, only : inrange 
  use usefulbits, only : g2,ms,mc,mbmb,mmu,mtau,small
  implicit none
  integer, intent(in) :: jj
  !--------------------------------------internal
  integer :: i
  double precision :: Mhi,GammaRat
  !---------------------------------------------

  do i=1,np(Hneut)
     
   Mhi=theo(jj)%particle(Hneut)%M(i)    
   if(theo(jj)%particle(Hneut)%Mc(i).ge.small) Mhi=theo(jj)%particle(Hneut)%Mc(i)
     
    theo(jj)%BR_hjss(i)    = 0.0D0
    theo(jj)%BR_hjcc(i)    = 0.0D0       
    theo(jj)%BR_hjbb(i)    = 0.0D0
    theo(jj)%BR_hjmumu(i)  = 0.0D0 
    theo(jj)%BR_hjtautau(i)= 0.0D0 
    theo(jj)%BR_hjWW(i)    = 0.0D0 
    theo(jj)%BR_hjZZ(i)    = 0.0D0 
    theo(jj)%BR_hjZga(i)   = 0.0D0
    theo(jj)%BR_hjgaga(i)  = 0.0D0
    theo(jj)%BR_hjgg(i)    = 0.0D0

    if( inrange(Mhi,'SMBR') )then
    
     GammaRat=theo(jj)%particle(Hneut)%GammaTot(i)/BRSM_GammaTot(Mhi)  
  
     if(theo(jj)%particle(Hneut)%GammaTot(i).gt.0.0D0)then
      theo(jj)%BR_hjss(i)    = ( g2(jj)%hjss_s(i)    +g2(jj)%hjss_p(i)    *invbsq(ms,  Mhi) ) *BRSM_Hss(Mhi)     /GammaRat
      theo(jj)%BR_hjcc(i)    = ( g2(jj)%hjcc_s(i)    +g2(jj)%hjcc_p(i)    *invbsq(mc,  Mhi) ) *BRSM_Hcc(Mhi)     /GammaRat
      theo(jj)%BR_hjbb(i)    = ( g2(jj)%hjbb_s(i)    +g2(jj)%hjbb_p(i)    *invbsq(mbmb,Mhi) ) *BRSM_Hbb(Mhi)     /GammaRat
      theo(jj)%BR_hjmumu(i)  = ( g2(jj)%hjmumu_s(i)  +g2(jj)%hjmumu_p(i)  *invbsq(mmu, Mhi) ) *BRSM_Hmumu(Mhi)   /GammaRat
      theo(jj)%BR_hjtautau(i)= ( g2(jj)%hjtautau_s(i)+g2(jj)%hjtautau_p(i)*invbsq(mtau,Mhi) ) *BRSM_Htautau(Mhi) /GammaRat

      theo(jj)%BR_hjWW(i)    = g2(jj)%hjWW(i)    *BRSM_HWW(Mhi)    /GammaRat 
      theo(jj)%BR_hjZZ(i)    = g2(jj)%hjZZ(i)    *BRSM_HZZ(Mhi)    /GammaRat 
      theo(jj)%BR_hjZga(i)   = g2(jj)%hjZga(i)   *BRSM_HZga(Mhi)   /GammaRat 
      theo(jj)%BR_hjgaga(i)  = g2(jj)%hjgaga(i)  *BRSM_Hgaga(Mhi)  /GammaRat  
      theo(jj)%BR_hjgg(i)    = g2(jj)%hjgg(i)    *BRSM_Hgg(Mhi)    /GammaRat  
     else
      write(*,*)'at jj=',jj,'i=',i
      write(*,*)'total decay width is less than or equal to zero:',theo(jj)%particle(Hneut)%GammaTot(i)
     endif     
    endif                          
   enddo 

 end subroutine br_from_g2_for_datapoint
 !*****************************************************************
 function invbsq(mf,mh)
  implicit none
  double precision,intent(in) :: mf,mh
  double precision :: invbsq
  if(mh>2.0D0*mf)then 
   invbsq=1.0D0/(1.0D0-4.0D0*(mf/mh)**2.0D0)
  else
   invbsq=0.0D0
  endif
 end function invbsq
 !*****************************************************************
 subroutine check_dataset
 ! checks each parameter point to determine whether the Higgs masses
 ! and branching ratios make sense
 ! Sets theo(jj)%gooddataset accordingly
 !*****************************************************************
  use usefulbits, only : theo,ndat,debug,np
  implicit none
  !--------------------------------------internal
  integer :: jj,x
  double precision :: testsumBR,testsumBR_t
  double precision,allocatable :: testBR(:)
  double precision :: fuzziness
  !---------------------------------------------
  fuzziness = 0.001D0
  !fuzziness = 100.0D0 ; write(*,*)'WARNING: fuzziness factor is far too high'

  testsumBR =0.0D0
  testsumBR_t  =0.0D0

  if(np(Hneut)>0)then
   allocate(testBR(np(Hneut)))

   ! testing to see if the dataset is ok
   do jj=1,ndat   
    testBR     =         theo(jj)%BR_hjss              &
            &      +     theo(jj)%BR_hjcc              &
            &      +     theo(jj)%BR_hjbb              &
            &      +     theo(jj)%BR_hjmumu            &
            &      +     theo(jj)%BR_hjtautau          &
            &      +     theo(jj)%BR_hjWW              &
            &      +     theo(jj)%BR_hjZZ              &
            &      +     theo(jj)%BR_hjZga             &
            &      +     theo(jj)%BR_hjgaga            &
            &      + sum(theo(jj)%BR_hjhihi,dim=2)
    testsumBR  = maxval(  testBR  ) 
    if( testsumBR   .gt.   1.0D0+fuzziness   )then     
     if(debug)write(*,*) 'warning: sum of BR for '//trim(adjustl(pdesc(Hneut)%long))//' at line number=',jj,'is',testsumBR
    endif
   enddo
   deallocate(testBR)
  endif 

  if(np(Hplus)>0)then
   allocate(testBR(np(Hplus)))

   do jj=1,ndat 
    testBR      =        theo(jj)%BR_Hpjcs             &
            &      +     theo(jj)%BR_Hpjcb             &
            &      +     theo(jj)%BR_Hpjtaunu          
    testsumBR = maxval(  testBR ) 


    testsumBR_t =        theo(jj)%BR_tWpb              &
            &      + sum(theo(jj)%BR_tHpjb,dim=1)

    if( testsumBR  .gt.   1.0D0+fuzziness   )then     
     if(debug)write(*,*) 'warning: sum of BR for '//trim(adjustl(pdesc(Hplus)%long))//' at line number=',jj,'is',testsumBR
    elseif( testsumBR_t    .gt.   1.0D0+fuzziness   )then     
     if(debug)write(*,*) 'warning: sum of BR for the top quark at jj=',jj,'is',testsumBR_t  
    endif 

   enddo
   deallocate(testBR)
  endif

  do jj=1,ndat
     theo(jj)%gooddataset=.True.      
  enddo 

  do x=1,ubound(np,dim=1)
   if(np(x)>0)then
    do jj=1,ndat 
     if(     minval(theo(jj)%particle(x)%M).lt.0.0D0)then
      theo(jj)%gooddataset=.False.
      if(debug)write(*,*) 'warning: negative mass for '//trim(adjustl(pdesc(x)%long))//' at line number=',jj,theo(jj)%particle(x)%M

    !elseif( testsumBR_hj  .gt.   (1.0D0+fuzziness)   )then !i.e. branching ratios for one of the Higgs add up to more than 1+fuzziness
    ! !theo(jj)%gooddataset=.False. 
     elseif( .not. (sum(theo(jj)%particle(x)%M).ge.0.0D0) )then   
      theo(jj)%gooddataset=.False.
      write(*,*) 'warning: mass is NaN for '//trim(adjustl(pdesc(x)%long))//' at line number=',jj,theo(jj)%particle(x)%M

     elseif(     minval(theo(jj)%particle(x)%GammaTot).lt.0.0D0)then
      theo(jj)%gooddataset=.False.
      if(debug)write(*,*) 'warning: negative total decay width for '//trim(adjustl(pdesc(x)%long))// &
                        & ' at line number=',jj,theo(jj)%particle(x)%GammaTot

    !elseif( testsumBR_hj  .gt.   (1.0D0+fuzziness)   )then !i.e. branching ratios for one of the Higgs add up to more than 1+fuzziness
    ! !theo(jj)%gooddataset=.False. 
     elseif( .not. (sum(theo(jj)%particle(x)%GammaTot).ge.0.0D0) )then   
      theo(jj)%gooddataset=.False.
      if(debug)write(*,*) 'warning: total decay width is NaN for '//trim(adjustl(pdesc(x)%long))// &
                        & ' at line number=',jj,theo(jj)%particle(x)%GammaTot

     endif 

    enddo
   endif
  enddo

 end subroutine check_dataset

 !*****************************************************************
 subroutine fill_theo_SM
 ! fills the Standard Model part of theo
 ! We do this here to save computational time - these  quantities will be
 ! needed a few times in subroutine calcfact_t1, so don't want to calculate them each time
 !************************************************************
  use theory_BRfunctions
  use theory_XS_SM_functions
  use usefulbits, only : ndat
  use S95tables, only : inrange
  implicit none
  !--------------------------------------internal
  integer :: n
  !----------------------------------------------
  
  if(np(Hneut)<1)stop'error in subroutine fill_theo_SM (np(Hneut))'

   do n=1,ndat        
    call fill_theo_SM_for_datapoint(n)                     
   enddo        
     
 end subroutine fill_theo_SM
 !*****************************************************************
 subroutine fill_theo_SM_for_datapoint(n)
 ! fills the Standard Model part of theo
 ! We do this here to save computational time - these  quantities will be
 ! needed a few times in subroutine calcfact_t1, so don't want to calculate them each time
 !************************************************************
  use theory_BRfunctions
  use theory_XS_SM_functions
  use usefulbits, only : theo,small
  use S95tables, only : inrange
  implicit none
  integer, intent(in) :: n
  !--------------------------------------internal
  integer :: i
  double precision :: Mhi
  !----------------------------------------------
         
  if(theo(n)%gooddataset) then             
   do i=1,np(Hneut)
  
   Mhi=theo(n)%particle(Hneut)%M(i)    
   if(theo(n)%particle(Hneut)%Mc(i).ge.small) Mhi=theo(n)%particle(Hneut)%Mc(i)
    
   if(inrange(Mhi,'SMBR'))then
    theo(n)%BR_HWW_SM(i)    = BRSM_HWW(Mhi) 
    theo(n)%BR_HZZ_SM(i)    = BRSM_HZZ(Mhi)    
    theo(n)%BR_Hbb_SM(i)    = BRSM_Hbb(Mhi)
    theo(n)%BR_Hcc_SM(i)    = BRSM_Hcc(Mhi)    
    theo(n)%BR_Hss_SM(i)    = BRSM_Hss(Mhi)        
    theo(n)%BR_Hmumu_SM(i)  = BRSM_Hmumu(Mhi)
    theo(n)%BR_Htautau_SM(i)= BRSM_Htautau(Mhi)
    theo(n)%BR_HZga_SM(i)   = BRSM_HZga(Mhi)
    theo(n)%BR_Hgaga_SM(i)  = BRSM_Hgaga(Mhi)
    theo(n)%BR_Hgg_SM(i)  = BRSM_Hgg(Mhi)    
    theo(n)%BR_Hjets_SM(i)  = BRSM_Hss(Mhi)+BRSM_Hcc(Mhi)+BRSM_Hbb(Mhi)+BRSM_Hgg(Mhi)
    theo(n)%GammaTot_SM(i)  = BRSM_GammaTot(Mhi)       
   else
    theo(n)%BR_HWW_SM(i)     = 0.0D0
    theo(n)%BR_HZZ_SM(i)     = 0.0D0
    theo(n)%BR_Hbb_SM(i)     = 0.0D0
    theo(n)%BR_Hcc_SM(i)     = 0.0D0
    theo(n)%BR_Hss_SM(i)     = 0.0D0    
    theo(n)%BR_Hmumu_SM(i)   = 0.0D0
    theo(n)%BR_Htautau_SM(i) = 0.0D0
    theo(n)%BR_HZga_SM(i)    = 0.0D0
    theo(n)%BR_Hgaga_SM(i)   = 0.0D0
    theo(n)%BR_Hgg_SM(i)     = 0.0D0    
    theo(n)%BR_Hjets_SM(i)   = 0.0D0
    theo(n)%GammaTot_SM(i)   = 0.0D0 
   endif

   if(inrange(Mhi,'TEV '))then             
    theo(n)%tev%XS_HZ_SM(i) = XS_tev_HZ_SM(Mhi)   
    theo(n)%tev%XS_HW_SM(i) = XS_tev_HW_SM(Mhi)
    theo(n)%tev%XS_H_SM(i)  = XS_tev_gg_H_SM(Mhi)+XS_tev_bb_H_SM(Mhi)  
    theo(n)%tev%XS_vbf_SM(i)= XS_tev_vbf_SM(Mhi)
    theo(n)%tev%XS_ttH_SM(i)= XS_tev_ttH_SM(Mhi)  
           
    theo(n)%tev%XS_Hb_SM(i)    = XS_tev_bg_Hb_SM(Mhi) 
    theo(n)%tev%XS_Hb_c1_SM(i) = XS_tev_bg_Hb_c1_SM(Mhi) 
    theo(n)%tev%XS_Hb_c2_SM(i) = XS_tev_bg_Hb_c2_SM(Mhi)    
    theo(n)%tev%XS_Hb_c3_SM(i) = XS_tev_bg_Hb_c3_SM(Mhi)  
    theo(n)%tev%XS_Hb_c4_SM(i) = XS_tev_bg_Hb_c4_SM(Mhi) 
   else
    theo(n)%tev%XS_HW_SM(i) = 0.0D0
    theo(n)%tev%XS_H_SM(i)  = 0.0D0
    theo(n)%tev%XS_HZ_SM(i) = 0.0D0
    theo(n)%tev%XS_vbf_SM(i)= 0.0D0
    theo(n)%tev%XS_ttH_SM(i)= 0.0D0
     
    theo(n)%tev%XS_Hb_SM(i)    = 0.0D0
    theo(n)%tev%XS_Hb_c1_SM(i) = 0.0D0
    theo(n)%tev%XS_Hb_c2_SM(i) = 0.0D0
    theo(n)%tev%XS_Hb_c3_SM(i) = 0.0D0 
    theo(n)%tev%XS_Hb_c4_SM(i) = 0.0D0 
   endif

   if(inrange(Mhi,'LHC7'))then
    theo(n)%lhc7%XS_HW_SM(i) = XS_lhc7_HW_SM(Mhi) 
    theo(n)%lhc7%XS_H_SM(i)  = XS_lhc7_gg_H_SM(Mhi) + XS_lhc7_bb_H_SM(Mhi) 
    theo(n)%lhc7%XS_HZ_SM(i) = XS_lhc7_HZ_SM(Mhi) 
    theo(n)%lhc7%XS_vbf_SM(i)= XS_lhc7_vbf_SM(Mhi) 
    theo(n)%lhc7%XS_ttH_SM(i)= XS_lhc7_ttH_SM(Mhi)
   else     
    theo(n)%lhc7%XS_HW_SM(i) = 0.0D0
    theo(n)%lhc7%XS_H_SM(i)  = 0.0D0
    theo(n)%lhc7%XS_HZ_SM(i) = 0.0D0
    theo(n)%lhc7%XS_vbf_SM(i)= 0.0D0             
    theo(n)%lhc7%XS_ttH_SM(i)= 0.0D0
   endif     

   if(inrange(Mhi,'LHC8'))then
    theo(n)%lhc8%XS_HW_SM(i) = XS_lhc8_HW_SM(Mhi) 
    theo(n)%lhc8%XS_H_SM(i)  = XS_lhc8_gg_H_SM(Mhi) + XS_lhc8_bb_H_SM(Mhi) 
    theo(n)%lhc8%XS_HZ_SM(i) = XS_lhc8_HZ_SM(Mhi) 
    theo(n)%lhc8%XS_vbf_SM(i)= XS_lhc8_vbf_SM(Mhi) 
    theo(n)%lhc8%XS_ttH_SM(i)= XS_lhc8_ttH_SM(Mhi)
   else     
    theo(n)%lhc8%XS_HW_SM(i) = 0.0D0
    theo(n)%lhc8%XS_H_SM(i)  = 0.0D0
    theo(n)%lhc8%XS_HZ_SM(i) = 0.0D0
    theo(n)%lhc8%XS_vbf_SM(i)= 0.0D0             
    theo(n)%lhc8%XS_ttH_SM(i)= 0.0D0
   endif     


  enddo
 endif                     
     
 end subroutine fill_theo_SM_for_datapoint
 !************************************************************
 subroutine XS_from_partR
 ! turn partonic cross section ratios in to hadronic cross section
 ! ratios
 ! Subroutine is complicated by the fact that if e.g. all
 ! the partR(n)%qq_hjW partonic cross section ratios are equal,
 ! just want to use this value for the hadronic cross section
 ! ratio and not lose any accuracy by combining with the tevS
 !************************************************************
  use usefulbits, only : ndat               
  use S95tables, only : inrange

  implicit none
  !--------------------------------------internal
  integer :: n
  !----------------------------------------------

  if(np(Hneut)<1)stop'error in subroutine XS_from_partR (np(Hneut))'

  do n=1,ndat        
   call XS_from_partR_for_datapoint(n)
  enddo            

 end subroutine XS_from_partR
 !******************************************************************
 subroutine XS_from_partR_for_datapoint(n)
 ! turn partonic cross section ratios in to hadronic cross section
 ! ratios
 ! Subroutine is complicated by the fact that if e.g. all
 ! the partR(n)%qq_hjW partonic cross section ratios are equal,
 ! just want to use this value for the hadronic cross section
 ! ratio and not lose any accuracy by combining with the tevS
 !************************************************************
  use usefulbits, only : allocate_hadroncolliderextras_parts, &
  &                      deallocate_hadroncolliderextras_parts
  use S95tables, only : inrange
  implicit none
  integer, intent(in) :: n
  !--------------------------------------internal
  integer :: i
  double precision :: Mhi
  logical :: simple_partR
  !----------------------------------------------

  call allocate_hadroncolliderextras_parts(tevS) 
  call allocate_hadroncolliderextras_parts(lhc7S) 
  call allocate_hadroncolliderextras_parts(lhc8S) 

  if(theo(n)%gooddataset) then             
   do i=1,np(Hneut)
    Mhi=theo(n)%particle(Hneut)%M(i)

    call fill_tevS(i,Mhi) 
    call fill_lhc7S(i,Mhi)     
    call fill_lhc8S(i,Mhi)     
	
	!this if loop is here to make sure partR(n)%qq_hjWp(1,i).eq.0.0D0 is taken care of
    if(partR(n)%qq_hjWp(1,i).eq.0.0D0)then
     simple_partR=.False.
    elseif( (( sum(abs( partR(n)%qq_hjWp(:,i) - partR(n)%qq_hjWp(1,i)))   &
&            + sum(abs( partR(n)%qq_hjWm(:,i) - partR(n)%qq_hjWp(1,i))) )/partR(n)%qq_hjWp(1,i)) .lt. 1.0D-5 )then
     simple_partR=.True.
    else
     simple_partR=.False.
    endif

   if(simple_partR)then
    theo(n)%tev%XS_hjW_ratio(i)=partR(n)%qq_hjWp(1,i)
    theo(n)%lhc7%XS_hjW_ratio(i)=partR(n)%qq_hjWp(1,i)
    theo(n)%lhc8%XS_hjW_ratio(i)=partR(n)%qq_hjWp(1,i)
   else
    theo(n)%tev%XS_hjW_ratio(i)=                          &
    &   sum( partR(n)%qq_hjWp(:,i)*tevS(1)%qq_hjWp(:,i) ) &
    & + sum( partR(n)%qq_hjWm(:,i)*tevS(1)%qq_hjWm(:,i) )
    theo(n)%lhc7%XS_hjW_ratio(i)=                          &
    &   sum( partR(n)%qq_hjWp(:,i)*lhc7S(1)%qq_hjWp(:,i) ) &
    & + sum( partR(n)%qq_hjWm(:,i)*lhc7S(1)%qq_hjWm(:,i) )
    theo(n)%lhc8%XS_hjW_ratio(i)=                          &
    &   sum( partR(n)%qq_hjWp(:,i)*lhc8S(1)%qq_hjWp(:,i) ) &
    & + sum( partR(n)%qq_hjWm(:,i)*lhc8S(1)%qq_hjWm(:,i) )
   endif

   theo(n)%tev%XS_hj_ratio(i)=                         &
    &        partR(n)%gg_hj(i)    *tevS(1)%gg_hj(i)     &
    & + sum( partR(n)%qq_hj(:,i)  *tevS(1)%qq_hj(:,i)   )
   theo(n)%lhc7%XS_hj_ratio(i)=                         &
    &        partR(n)%gg_hj(i)    *lhc7S(1)%gg_hj(i)     &
    & + sum( partR(n)%qq_hj(:,i)  *lhc7S(1)%qq_hj(:,i)   ) 
   theo(n)%lhc8%XS_hj_ratio(i)=                         &
    &        partR(n)%gg_hj(i)    *lhc8S(1)%gg_hj(i)     &
    & + sum( partR(n)%qq_hj(:,i)  *lhc8S(1)%qq_hj(:,i)   ) 

   !this if loop is here to make sure partR(n)%qq_hjZ(1,i).eq.0.0D0 is taken care of
   if(partR(n)%qq_hjZ(1,i).eq.0.0D0)then
    simple_partR=.False.
   elseif( (abs(sum( partR(n)%qq_hjZ(:,i) - partR(n)%qq_hjZ(1,i)))/partR(n)%qq_hjZ(1,i)).lt. 1.0D-5 )then
    simple_partR=.True.
   else
    simple_partR=.False.
   endif

   if(  simple_partR  )then 
    theo(n)%tev%XS_hjZ_ratio(i) =  partR(n)%qq_hjZ(1,i)
    if(partR(n)%gg_hjZ(i) .le.0.0D0)then
     theo(n)%lhc7%XS_hjZ_ratio(i)=  partR(n)%qq_hjZ(1,i)
     theo(n)%lhc8%XS_hjZ_ratio(i)=  partR(n)%qq_hjZ(1,i)
    else
     theo(n)%lhc7%XS_hjZ_ratio(i)=                     &
     &        partR(n)%gg_hjZ(i)   *lhc7S(1)%gg_hjZ(i)    &
     & + sum( partR(n)%qq_hjZ(:,i) *lhc7S(1)%qq_hjZ(:,i)  )
     theo(n)%lhc8%XS_hjZ_ratio(i)=                     &
     &        partR(n)%gg_hjZ(i)   *lhc8S(1)%gg_hjZ(i)    &
     & + sum( partR(n)%qq_hjZ(:,i) *lhc8S(1)%qq_hjZ(:,i)  )
    endif
   else
    theo(n)%tev%XS_hjZ_ratio(i)=                     &
    &        partR(n)%gg_hjZ(i)   *tevS(1)%gg_hjZ(i)    &
    & + sum( partR(n)%qq_hjZ(:,i) *tevS(1)%qq_hjZ(:,i)  )
    theo(n)%lhc7%XS_hjZ_ratio(i)=                     &
    &        partR(n)%gg_hjZ(i)   *lhc7S(1)%gg_hjZ(i)    &
    & + sum( partR(n)%qq_hjZ(:,i) *lhc7S(1)%qq_hjZ(:,i)  )
    theo(n)%lhc8%XS_hjZ_ratio(i)=                     &
    &        partR(n)%gg_hjZ(i)   *lhc8S(1)%gg_hjZ(i)    &
    & + sum( partR(n)%qq_hjZ(:,i) *lhc8S(1)%qq_hjZ(:,i)  )
   endif

   theo(n)%tev%XS_hjb_ratio(i) = partR(n)%bg_hjb(i)
   theo(n)%lhc7%XS_hjb_ratio(i)= partR(n)%bg_hjb(i) 
   theo(n)%lhc8%XS_hjb_ratio(i)= partR(n)%bg_hjb(i) 

   enddo
  endif                       

  call deallocate_hadroncolliderextras_parts(lhc8S)
  call deallocate_hadroncolliderextras_parts(lhc7S)
  call deallocate_hadroncolliderextras_parts(tevS)

 end subroutine XS_from_partR_for_datapoint
 !******************************************************************   
   
 subroutine fill_tevS(j,Mhj)
 ! fills the elements of tevS using the functions in module theory_colliderSfunctions
 !************************************************************
  use theory_colliderSfunctions
  use theory_XS_SM_functions
  use S95tables, only : inrange
  implicit none
  !--------------------------------------internal
  integer :: j
  double precision :: Mhj
  !----------------------------------------------

  if(inrange(Mhj,'TEV '))then
    tevS(1)%qq_hjWp(1,j)=tev_rHWpm_udb(Mhj)
    tevS(1)%qq_hjWp(2,j)=tev_rHWpm_csb(Mhj) 
   
    tevS(1)%qq_hjWm(1,j)=tev_rHWpm_dub(Mhj)
    tevS(1)%qq_hjWm(2,j)=tev_rHWpm_scb(Mhj)
  
    !We now have a new gg->H SM function: Should use XS functions instead of r's
    !For cross check with OB code changed this temporarily!
    tevS(1)%gg_hj(j)=tev_rH_gg(Mhj) 
    !tevS(1)%gg_hj(j)  =XS_tev_gg_H_SM(Mhj)/(XS_tev_gg_H_SM(Mhj)+XS_tev_bb_H_SM(Mhj)) 
    tevS(1)%qq_hj(:,j)=0.0D0
    tevS(1)%qq_hj(5,j)=tev_rH_bb(Mhj)
    !tevS(1)%qq_hj(5,j)=XS_tev_bb_H_SM(Mhj)/(XS_tev_gg_H_SM(Mhj)+XS_tev_bb_H_SM(Mhj)) 

    tevS(1)%gg_hjZ(j)=0.0D0
  
    tevS(1)%qq_hjZ(1,j)=tev_rHZ_ddb(Mhj)
    tevS(1)%qq_hjZ(2,j)=tev_rHZ_uub(Mhj)
    tevS(1)%qq_hjZ(3,j)=tev_rHZ_ssb(Mhj)
    tevS(1)%qq_hjZ(4,j)=tev_rHZ_ccb(Mhj)
    tevS(1)%qq_hjZ(5,j)=tev_rHZ_bbb(Mhj)
      

    if(abs(sum(tevS(1)%qq_hjWp(:,j))+sum(tevS(1)%qq_hjWm(:,j)) - 1.0D0) .gt. 1.0D-2)then 
     stop 'error in fill_tevS (a)'
    elseif(abs(tevS(1)%gg_hj(j)+sum(tevS(1)%qq_hj(:,j)) - 1.0D0) .gt. 1.0D-2)then
     stop 'error in fill_tevS (b)'
    elseif(abs(tevS(1)%gg_hjZ(j)+sum(tevS(1)%qq_hjZ(:,j)) - 1.0D0) .gt. 1.0D-2)then
     stop 'error in fill_tevS (c)'
    endif 
  else
    tevS(1)%qq_hjWp(:,j)=0.0D0
    tevS(1)%qq_hjWm(:,j)=0.0D0
    tevS(1)%gg_hj(j)=0.0D0
    tevS(1)%qq_hj(:,j)=0.0D0
    tevS(1)%gg_hjZ(j)=0.0D0
    tevS(1)%qq_hjZ(:,j)=0.0D0
  endif

 end subroutine fill_tevS 
 !******************************************************************  

 subroutine fill_lhc7S(j,Mhj)
 ! fills the elements of lhc7S using the functions in module theory_colliderSfunctions
 !************************************************************
  use theory_colliderSfunctions
  use theory_XS_SM_functions
  use usefulbits, only : vsmall
  use S95tables, only : inrange
  implicit none
  !--------------------------------------internal
  integer :: j
  double precision :: Mhj
  !----------------------------------------------
  
  if(inrange(Mhj,'LHC7'))then
    lhc7S(1)%gg_hj(j)=LHC7_rH_gg(Mhj)
    
    lhc7S(1)%qq_hj(:,j)=0.0D0
    lhc7S(1)%qq_hj(5,j)=LHC7_rH_bb(Mhj)

    if(abs(lhc7S(1)%gg_hj(j)+sum(lhc7S(1)%qq_hj(:,j)) - 1.0D0) .gt. 1.0D-2)then
     stop 'error in fill_lhc7S (b)'
    endif  

    if(XS_lhc7_HW_SM(Mhj).lt.vsmall)then
      lhc7S(1)%qq_hjWp(:,j)=0.0D0
      lhc7S(1)%qq_hjWm(:,j)=0.0D0
    else
      lhc7S(1)%qq_hjWp(1,j)=LHC7_rHWp_udb(Mhj)
      lhc7S(1)%qq_hjWp(2,j)=LHC7_rHWp_csb(Mhj) 
   
      lhc7S(1)%qq_hjWm(1,j)=LHC7_rHWm_dub(Mhj)
      lhc7S(1)%qq_hjWm(2,j)=LHC7_rHWm_scb(Mhj)
  
      if(abs(sum(lhc7S(1)%qq_hjWp(:,j))+sum(lhc7S(1)%qq_hjWm(:,j)) - 1.0D0) .gt. 1.0D-2)then 
       stop 'error in fill_lhc7S (a)'
      endif 
    endif
    if(XS_lhc7_HZ_SM(Mhj).lt.vsmall)then
      lhc7S(1)%gg_hjZ(j)=0.0D0
      lhc7S(1)%qq_hjZ(:,j)=0.0D0
    else
      lhc7S(1)%gg_hjZ(j)=LHC7_rHZ_gg(Mhj)
  
      lhc7S(1)%qq_hjZ(1,j)=LHC7_rHZ_ddb(Mhj)
      lhc7S(1)%qq_hjZ(2,j)=LHC7_rHZ_uub(Mhj)
      lhc7S(1)%qq_hjZ(3,j)=LHC7_rHZ_ssb(Mhj)
      lhc7S(1)%qq_hjZ(4,j)=LHC7_rHZ_ccb(Mhj)
      lhc7S(1)%qq_hjZ(5,j)=LHC7_rHZ_bbb(Mhj)

      if(abs(lhc7S(1)%gg_hjZ(j)+sum(lhc7S(1)%qq_hjZ(:,j)) - 1.0D0) .gt. 1.0D-2)then
       stop 'error in fill_lhc7S (c)'
      endif 
    endif   
  else
    lhc7S(1)%qq_hjWp(:,j)=0.0D0
    lhc7S(1)%qq_hjWm(:,j)=0.0D0
    lhc7S(1)%gg_hj(j)=0.0D0
    lhc7S(1)%qq_hj(:,j)=0.0D0
    lhc7S(1)%gg_hjZ(j)=0.0D0
    lhc7S(1)%qq_hjZ(:,j)=0.0D0
  endif

 end subroutine fill_lhc7S 
 !******************************************************************  

 subroutine fill_lhc8S(j,Mhj)
 ! fills the elements of lhc8S using the functions in module theory_colliderSfunctions
 !************************************************************
  use theory_colliderSfunctions
  use theory_XS_SM_functions
  use usefulbits, only : vsmall
  use S95tables, only : inrange
  implicit none
  !--------------------------------------internal
  integer :: j
  double precision :: Mhj
  !----------------------------------------------
  
  if(inrange(Mhj,'LHC8'))then
    lhc8S(1)%gg_hj(j)=LHC8_rH_gg(Mhj)
    
    lhc8S(1)%qq_hj(:,j)=0.0D0
    lhc8S(1)%qq_hj(5,j)=LHC8_rH_bb(Mhj)

    if(abs(lhc8S(1)%gg_hj(j)+sum(lhc8S(1)%qq_hj(:,j)) - 1.0D0) .gt. 1.0D-2)then
     stop 'error in fill_lhc8S (b)'
    endif  

    if(XS_lhc8_HW_SM(Mhj).lt.vsmall)then
      lhc8S(1)%qq_hjWp(:,j)=0.0D0
      lhc8S(1)%qq_hjWm(:,j)=0.0D0
    else
      lhc8S(1)%qq_hjWp(1,j)=LHC8_rHWp_udb(Mhj)
      lhc8S(1)%qq_hjWp(2,j)=LHC8_rHWp_csb(Mhj) 
   
      lhc8S(1)%qq_hjWm(1,j)=LHC8_rHWm_dub(Mhj)
      lhc8S(1)%qq_hjWm(2,j)=LHC8_rHWm_scb(Mhj)
  
      if(abs(sum(lhc8S(1)%qq_hjWp(:,j))+sum(lhc8S(1)%qq_hjWm(:,j)) - 1.0D0) .gt. 1.0D-2)then 
       stop 'error in fill_lhc8S (a)'
      endif 
    endif
    if(XS_lhc8_HZ_SM(Mhj).lt.vsmall)then
      lhc8S(1)%gg_hjZ(j)=0.0D0
      lhc8S(1)%qq_hjZ(:,j)=0.0D0
    else
      lhc8S(1)%gg_hjZ(j)=LHC8_rHZ_gg(Mhj)
  
      lhc8S(1)%qq_hjZ(1,j)=LHC8_rHZ_ddb(Mhj)
      lhc8S(1)%qq_hjZ(2,j)=LHC8_rHZ_uub(Mhj)
      lhc8S(1)%qq_hjZ(3,j)=LHC8_rHZ_ssb(Mhj)
      lhc8S(1)%qq_hjZ(4,j)=LHC8_rHZ_ccb(Mhj)
      lhc8S(1)%qq_hjZ(5,j)=LHC8_rHZ_bbb(Mhj)

      if(abs(lhc8S(1)%gg_hjZ(j)+sum(lhc8S(1)%qq_hjZ(:,j)) - 1.0D0) .gt. 1.0D-2)then
       stop 'error in fill_lhc8S (c)'
      endif 
    endif   
  else
    lhc8S(1)%qq_hjWp(:,j)=0.0D0
    lhc8S(1)%qq_hjWm(:,j)=0.0D0
    lhc8S(1)%gg_hj(j)=0.0D0
    lhc8S(1)%qq_hj(:,j)=0.0D0
    lhc8S(1)%gg_hjZ(j)=0.0D0
    lhc8S(1)%qq_hjZ(:,j)=0.0D0
  endif

 end subroutine fill_lhc8S 
 !******************************************************************  


end module theo_manip
!******************************************************************
