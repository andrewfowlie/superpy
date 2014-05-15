! This file is part of HiggsBounds
!  -KW
!******************************************************************
module S95tables
!******************************************************************      
 use S95tables_type1      
 use S95tables_type2  
 use usefulbits, only: Hneut,Hplus
 implicit none
 private 
 public ::  &
         & S95_t1,S95_t2, &
         & setup_S95tables,deallocate_S95tables, &
         & calcfact,outputproc, &
         & check_against_bound, &
         & convolve_chisq_with_gaussian, &
         & S95_t1_or_S95_t2_idfromelementnumber, &
         & f_from_t2,f_from_t3,f_from_slices_t2, &
         & WhichColliderString,get_collider_element_number, &
         & inrange, deallocate_Exptranges
  

 integer :: ntable1,ntable2

 !-----------------------------------------------------------------------------------
 ! delta_M*_* determine how close in mass particles should be before their masses are combined

 ! nb. delta*_* is only used for tables of type 1 at the moment

  ! for some analyses, S95_t1(x)%deltax has already been set in S95tables_type1.f90
  ! for analyses where S95_t1(x)%deltax is *not* has already set:

  ! we set 
  !                       LEP neutral Higgs tables to have deltax=delta_Mh_LEP
  !                       Tevatron neutral Higgs tables to have deltax=delta_Mh_TEV
  !                       LHC neutral Higgs tables to have deltax=delta_Mh_LHC
  !                       LEP charged Higgs tables to have deltax=delta_Mhplus_LEP
  !                       Tevatron charged Higgs tables to have deltax=delta_Mhplus_TEV
  !                       LHC charged Higgs tables to have deltax=delta_Mhplus_LHC
  !                       other LEP tables to have deltax=delta_M_LEP_default
  !                       other Tevatron tables to have deltax=delta_M_TEV_default
  !                       other LHC tables to have deltax=delta_M_LHC_default
  ! setting delta_M*_LHC and/or delta_M*_TEV and/or delta_M*_LEP to zero turns this feature off
  ! DO NOT CHANGE THESE VALUES BEFORE READING THE MANUAL
  ! even when it is appropriate to add the cross sections, we would recommend using 
  !    delta_Mh_LEP<=2.0D0, delta_Mh_TEV<=10.0, delta_Mh_LHC<=10.0
 double precision, parameter :: delta_Mh_LEP=0.0D0  
 double precision, parameter :: delta_Mh_TEV=10.0D0 
 double precision, parameter :: delta_Mh_LHC=10.0D0 
 double precision, parameter :: delta_Mhplus_LEP=0.0D0  
 double precision, parameter :: delta_Mhplus_TEV=0.0D0 
 double precision, parameter :: delta_Mhplus_LHC=0.0D0 
 double precision, parameter :: delta_M_LEP_default=0.0D0 !where delta_M_LEP/TEV/LHC is not specified 
 double precision, parameter :: delta_M_TEV_default=0.0D0 !these values will be used
 double precision, parameter :: delta_M_LHC_default=0.0D0 !

 !double precision, parameter :: delta_Mh_LEP=15.0D0 !crazy values - for debugging only 
 !double precision, parameter :: delta_Mh_TEV=15.0D0 !crazy(ish) values - for debugging only
 !double precision, parameter :: delta_Mh_LHC=15.0D0 !crazy(ish) values - for debugging only

 !-----------------------------------------------------------------------------------
 ! Use the SM expected channel contributions to weight the channels in the likeness test
 logical :: use_weight = .True.
 ! eps determines how strict the Standard Model test is 
 double precision, parameter :: eps=2.0D-2
 !double precision, parameter :: eps=1.0D3 !crazy value - for debugging only

 !-----------------------------------------------------------------------------------

 !table type 1-----------------------------      
 type(table1),allocatable :: S95_t1(:)
 !------------------------------------------
 
 !table type 2------------------------------
 type(table2),allocatable :: S95_t2(:)                              
 !------------------------------------------  
  character(LEN=4),parameter :: colliders(4) = (/'LEP ','TEV ','LHC7','LHC8'/)
  !-------------------------------------------  
  double precision, allocatable :: Exptrange_Mhmin_forSMdecays(:), Exptrange_Mhmax_forSMdecays(:)                         
  double precision, allocatable :: Exptrange_Mhmin_forSMXS(:),     Exptrange_Mhmax_forSMXS(:)

 contains

 !**********************************************************
 subroutine setup_S95tables
 ! Allocates  and calls subroutines to fill S95_t1, S95_t2 
 ! (which store the experimental data)
 ! Sets delta_M_TEV,delta_M_LEP,delta_M_LHC which govern how close Higgs
 ! need to be in mass before HiggsBounds combines their cross sections
 !**********************************************************
  use usefulbits, only : debug,np,not_a_particle
  use theory_BRfunctions, only : BRSMt1Mhmax,BRSMt1Mhmin 
  use theory_XS_SM_functions, only : tevXS_SM_functions_xmin, tevXS_SM_functions_xmax, &
                                   & lhc7XS_SM_functions_xmin,lhc7XS_SM_functions_xmax, &
                                   & lhc8XS_SM_functions_xmin,lhc8XS_SM_functions_xmax
                                   
  implicit none      
  !-----------------------------------internal 
  integer :: i,c
  double precision, allocatable :: max_expt_delta_Mh(:)
  double precision, allocatable :: Expttables_Mhmin_forSMXS(:),Expttables_Mhmax_forSMXS(:)
  double precision, allocatable :: Expttables_Mhmin_forSMdecays(:),Expttables_Mhmax_forSMdecays(:)
  double precision, allocatable :: delta_x_default(:,:) 

  ! these numbers have to be changed appropriately every time a table is added
  ! or taken away:  
  ntable1=133
  ntable2=17                ! table type 2 involves 2 variables
  
  allocate(S95_t1(ntable1))            
  allocate(S95_t2(ntable2))      

  call initializetables_type1_blank(S95_t1)
  call initializetables_type2_blank(S95_t2)

  call initializetables1(S95_t1)
  call initializetables2(S95_t2)    

  
  allocate( Expttables_Mhmin_forSMXS(  size(colliders,dim=1) ) )
  allocate( Expttables_Mhmax_forSMXS(  size(colliders,dim=1) ) )
  allocate( Exptrange_Mhmin_forSMXS(   size(colliders,dim=1) ) ) 
  allocate( Exptrange_Mhmax_forSMXS(   size(colliders,dim=1) ) )
  allocate( Expttables_Mhmin_forSMdecays(  size(colliders,dim=1) ) )
  allocate( Expttables_Mhmax_forSMdecays(  size(colliders,dim=1) ) )
  allocate( Exptrange_Mhmin_forSMdecays(   size(colliders,dim=1) ) ) 
  allocate( Exptrange_Mhmax_forSMdecays(   size(colliders,dim=1) ) )
  allocate( max_expt_delta_Mh( size(colliders,dim=1) ) )

  ! will pick up on any typos in S95_t1%expt
  do i=lbound(S95_t1,dim=1),ubound(S95_t1,dim=1)
   if(WhichColliderElement(S95_t1(i)%expt, S95_t1(i)%energy).eq.0)then
     write(*,*)'~'//trim(adjustl(S95_t1(i)%expt))//'~ is not a valid experiment name'
     stop 'error in setup_S95tables (a)'
   endif
  enddo

  ! will pick up on any typos in S95_t2%expt
  do i=lbound(S95_t2,dim=1),ubound(S95_t2,dim=1)
   if(WhichColliderElement(S95_t2(i)%expt, S95_t2(i)%energy).eq.0)then
     write(*,*)'~'//trim(adjustl(S95_t2(i)%expt))//'~ is not a valid experiment name'
     stop 'error in setup_S95tables (b)'
   endif
  enddo 

  ! checks that none of the id's are repeated in S95_t1
  do i=lbound(S95_t1,dim=1),ubound(S95_t1,dim=1)
   if( count(S95_t1%id.eq.S95_t1(i)%id) &
      +count(S95_t2%id.eq.S95_t1(i)%id).ne.1)then
    write(*,*)'the id',S95_t1(i)%id,'is repeated'
    stop 'error in setup_S95tables (c1)'
   endif
  enddo 
  
  ! checks that none of the id's are repeated in S95_t2
  do i=lbound(S95_t2,dim=1),ubound(S95_t2,dim=1)
   if( count(S95_t1%id.eq.S95_t2(i)%id)&
      +count(S95_t2%id.eq.S95_t2(i)%id).ne.1)then
    write(*,*)'the id',S95_t2(i)%id,'is repeated'
    stop 'error in setup_S95tables (c2)'
   endif
  enddo 

  !check to make sure that S95_t1(i)%particle_x are all particles
  do i=lbound(S95_t1,dim=1),ubound(S95_t1,dim=1)
   if(S95_t1(i)%particle_x.eq.not_a_particle)then
    write(*,*)S95_t1(i)%id,'particle_x=not_a_particle.'
    stop 'error in setup_S95tables (d1)'  
   endif
  enddo

  !check to make sure that S95_t2(i)%particle_x2 are all particles
  do i=lbound(S95_t2,dim=1),ubound(S95_t2,dim=1)
   if(S95_t2(i)%particle_x2.eq.not_a_particle)then
    write(*,*)S95_t2(i)%id,'particle_x2=not_a_particle.'
    stop 'error in setup_S95tables (d2)'  
   endif
  enddo
  
  ! looks for the min and max values of Mh (neutral Higgs) in each type of tables
  ! will be used in input.f90 to work out which SM para need to be calculated  
  
  ! initial (impossible) values
  Expttables_Mhmin_forSMXS=1.0D6
  Expttables_Mhmax_forSMXS=0.0D0

  do i=lbound(S95_t1,dim=1),ubound(S95_t1,dim=1)
    if(S95_t1(i)%particle_x.eq.Hneut)then
     do c=1,ubound(colliders,dim=1)
      if(WhichColliderElement(S95_t1(i)%expt, S95_t1(i)%energy).eq.c)then
       if(S95_t1(i)%xmax.gt.Expttables_Mhmax_forSMXS(c))Expttables_Mhmax_forSMXS(c)=S95_t1(i)%xmax
       if(S95_t1(i)%xmin.lt.Expttables_Mhmin_forSMXS(c))Expttables_Mhmin_forSMXS(c)=S95_t1(i)%xmin 
      endif 
     enddo
    endif
  enddo  

  do i=lbound(S95_t2,dim=1),ubound(S95_t2,dim=1)
   if(S95_t2(i)%particle_x2.eq.Hneut)then
    do c=1,ubound(colliders,dim=1)
     if(WhichColliderElement(S95_t2(i)%expt, S95_t2(i)%energy).eq.c)then
       if(S95_t2(i)%xmax2.gt.Expttables_Mhmax_forSMdecays(c))Expttables_Mhmax_forSMdecays(c)=S95_t2(i)%xmax2
       if(S95_t2(i)%xmin2.lt.Expttables_Mhmin_forSMdecays(c))Expttables_Mhmin_forSMdecays(c)=S95_t2(i)%xmin2 
     endif
    enddo
   endif
  enddo 

  Expttables_Mhmin_forSMdecays=Expttables_Mhmin_forSMXS
  Expttables_Mhmax_forSMdecays=Expttables_Mhmax_forSMXS

  do i=lbound(S95_t2,dim=1),ubound(S95_t2,dim=1)
   if(S95_t2(i)%particle_x1.eq.Hneut)then
    do c=1,ubound(colliders,dim=1)
     if(WhichColliderElement(S95_t2(i)%expt, S95_t2(i)%energy).eq.c)then
       if(S95_t2(i)%xmax1.gt.Expttables_Mhmax_forSMdecays(c))Expttables_Mhmax_forSMdecays(c)=S95_t2(i)%xmax1 
       if(S95_t2(i)%xmin1.lt.Expttables_Mhmin_forSMdecays(c))Expttables_Mhmin_forSMdecays(c)=S95_t2(i)%xmin1  
       ! Needs_M2_gt_2M1 is true only for processes involving hj->hihi.
       ! Therefore, SM production cross sections will not be needed for hi.
       if(.not.S95_t2(i)%needs_M2_gt_2M1)then
         if(S95_t2(i)%xmax1.gt.Expttables_Mhmax_forSMXS(c))Expttables_Mhmax_forSMXS(c)=S95_t2(i)%xmax1
         if(S95_t2(i)%xmin1.lt.Expttables_Mhmin_forSMXS(c))Expttables_Mhmin_forSMXS(c)=S95_t2(i)%xmin1      
       else
         if(S95_t2(i)%xmax2.gt.Expttables_Mhmax_forSMXS(c))Expttables_Mhmax_forSMXS(c)=S95_t2(i)%xmax2
         if(S95_t2(i)%xmin2.lt.Expttables_Mhmin_forSMXS(c))Expttables_Mhmin_forSMXS(c)=S95_t2(i)%xmin2            
       endif 
     endif 
    enddo
   endif
  enddo   

  ! now we set delta_x

  if(delta_M_LEP_default.gt.2.1d0) write(*,*)'WARNING: delta_M_LEP_default.gt.2.1d0'
  if(delta_M_TEV_default.gt.10.1d0)write(*,*)'WARNING: delta_M_TEV_default.gt.10.1d0'
  if(delta_M_LHC_default.gt.10.1d0)write(*,*)'WARNING: delta_M_LHC_default.gt.10.1d0'
  if(delta_Mh_LEP       .gt.2.1d0) write(*,*)'WARNING: delta_Mh_LEP.gt.2.1d0'
  if(delta_Mh_TEV       .gt.10.1d0)write(*,*)'WARNING: delta_Mh_TEV.gt.10.1d0'
  if(delta_Mh_LHC       .gt.10.1d0)write(*,*)'WARNING: delta_Mh_LHC.gt.10.1d0'
  if(delta_Mhplus_LEP   .gt.2.1d0) write(*,*)'WARNING: delta_Mhplus_LEP.gt.2.1d0'
  if(delta_Mhplus_TEV   .gt.10.1d0)write(*,*)'WARNING: delta_Mhplus_TEV.gt.10.1d0'
  if(delta_Mhplus_LHC   .gt.10.1d0)write(*,*)'WARNING: delta_Mhplus_LHC.gt.10.1d0'

  allocate( delta_x_default(size(np,dim=1),size(colliders,dim=1)) )

  ! fill delta_x_default
  do c=1,ubound(colliders,dim=1)
   if(c.eq.get_collider_element_number('LEP'))then ! for some reason, gfortran didn't like having a case statement here
     delta_x_default(:,c)    =delta_M_LEP_default
     delta_x_default(Hneut,c)=delta_Mh_LEP
     delta_x_default(Hplus,c)=delta_Mhplus_LEP
   elseif(c.eq.get_collider_element_number('TEV'))then
     delta_x_default(:,c)    =delta_M_TEV_default
     delta_x_default(Hneut,c)=delta_Mh_TEV
     delta_x_default(Hplus,c)=delta_Mhplus_TEV
   elseif(c.eq.get_collider_element_number('LHC7'))then
     delta_x_default(:,c)    =delta_M_LHC_default
     delta_x_default(Hneut,c)=delta_Mh_LHC
     delta_x_default(Hplus,c)=delta_Mhplus_LHC
   elseif(c.eq.get_collider_element_number('LHC8'))then
     delta_x_default(:,c)    =delta_M_LHC_default
     delta_x_default(Hneut,c)=delta_Mh_LHC
     delta_x_default(Hplus,c)=delta_Mhplus_LHC
   else
     stop'error in subroutine setup_S95tables'
   endif
  enddo

  do i=lbound(S95_t1,dim=1),ubound(S95_t1,dim=1)
   if(S95_t1(i)%deltax.lt.-0.5D0)then !i.e. deltax has not been set yet
    S95_t1(i)%deltax = delta_x_default(S95_t1(i)%particle_x,WhichColliderElement(S95_t1(i)%expt, S95_t1(i)%energy)) 
   endif
  enddo
   
  do i=lbound(S95_t2,dim=1),ubound(S95_t2,dim=1)
   if(S95_t2(i)%deltax.lt.-0.5D0)then !i.e. deltax has not been set yet
    S95_t2(i)%deltax = delta_x_default(S95_t2(i)%particle_x2,WhichColliderElement(S95_t2(i)%expt, S95_t2(i)%energy)) 
   endif
  enddo  
  
  ! finds the maximum delta_Mh for the each set of tables
  ! will be used in theo_SM.f90 to work out which SM para need to be calculated
  max_expt_delta_Mh= -1.0D0
  do i=lbound(S95_t1,dim=1),ubound(S95_t1,dim=1)
   do c=1,ubound(colliders,dim=1)
    if(WhichColliderElement(S95_t1(i)%expt, S95_t1(i)%energy).eq.c)then
     if( S95_t1(i)%particle_x.eq.Hneut)then
      if(S95_t1(i)%deltax    .gt.max_expt_delta_Mh(c))then
       max_expt_delta_Mh(c)=S95_t1(i)%deltax
      endif 
    endif
   endif
   enddo
  enddo  
  do i=lbound(S95_t2,dim=1),ubound(S95_t2,dim=1)
   do c=1,ubound(colliders,dim=1)
    if(WhichColliderElement(S95_t2(i)%expt, S95_t2(i)%energy).eq.c)then
     if( S95_t2(i)%particle_x2.eq.Hneut)then !note, this means only cross sections for different particle_x2 will be combined
      if(S95_t2(i)%deltax    .gt.max_expt_delta_Mh(c))then
       max_expt_delta_Mh(c)=S95_t2(i)%deltax
      endif 
    endif
   endif
   enddo
  enddo

 
  if(debug)write(*,*)'max_expt_delta_Mh',max_expt_delta_Mh

  Exptrange_Mhmin_forSMXS = max(Expttables_Mhmin_forSMXS  - max_expt_delta_Mh,0.0D0)
  Exptrange_Mhmax_forSMXS =     Expttables_Mhmax_forSMXS  + max_expt_delta_Mh

  Exptrange_Mhmin_forSMdecays = max(Expttables_Mhmin_forSMdecays  - max_expt_delta_Mh,0.0D0)
  Exptrange_Mhmax_forSMdecays =     Expttables_Mhmax_forSMdecays  + max_expt_delta_Mh

  !we need tevXS_SM_functions to have a big enough range to cover the tables
  if(Exptrange_Mhmax_forSMXS(get_collider_element_number('TEV')).gt.tevXS_SM_functions_xmax)then
   stop'need to extend upper range of tevXS_SM_functions or reduce delta_M_TEV'
  endif

  if(Exptrange_Mhmin_forSMXS(get_collider_element_number('TEV')).lt.tevXS_SM_functions_xmin)then
   write(*,*)Exptrange_Mhmin_forSMXS(get_collider_element_number('TEV')),tevXS_SM_functions_xmin
   stop'need to extend lower range of tevXS_SM_functions'
  endif

  !we need lhc7XS_SM_functions to have a big enough range to cover the tables
  if(Exptrange_Mhmax_forSMXS(get_collider_element_number('LHC7')).gt.lhc7XS_SM_functions_xmax)then
   stop'need to extend upper range of lhc7XS_SM_functions or reduce delta_M_LHC'
  endif

  if(Exptrange_Mhmin_forSMXS(get_collider_element_number('LHC7')).lt.lhc7XS_SM_functions_xmin)then
   stop'need to extend lower range of lhc7XS_SM_functions'
  endif

  !we need lhc7XS_SM_functions to have a big enough range to cover the tables
  if(Exptrange_Mhmax_forSMXS(get_collider_element_number('LHC8')).gt.lhc8XS_SM_functions_xmax)then
   stop'need to extend upper range of lhc8XS_SM_functions or reduce delta_M_LHC'
  endif

  if(Exptrange_Mhmin_forSMXS(get_collider_element_number('LHC8')).lt.lhc8XS_SM_functions_xmin)then
   stop'need to extend lower range of lhc8XS_SM_functions'
  endif

  
  ! we need the branching ratios for all the colliders
  if(    maxval(Exptrange_Mhmax_forSMdecays).gt.BRSMt1Mhmax)then
   stop'need to extend upper range of BRfunctions or reduce delta_M_(LEP/TEV)'
  elseif(minval(Exptrange_Mhmin_forSMdecays).lt.BRSMt1Mhmin)then
   write(*,*)'hello',minval(Exptrange_Mhmin_forSMdecays),BRSMt1Mhmin
   stop'need to extend lower range of BRfunctions'
  endif

  deallocate(Expttables_Mhmin_forSMXS)
  deallocate(Expttables_Mhmax_forSMXS)
  deallocate(Expttables_Mhmin_forSMdecays)
  deallocate(Expttables_Mhmax_forSMdecays)

  deallocate(max_expt_delta_Mh)
  deallocate(delta_x_default)

 end subroutine setup_S95tables      
 !**********************************************************  
 function inrange(mass,str)
 ! mass is the neutral Higgs mass to be checked
 ! str indicates which range it should be checked against:
 !   str should be either 'SMBR' or one of the elements of the array 'colliders'
 !
 ! The function returns .True. if the mass is in the appropriate range.
 ! If str='SMBR', this range is 
 !   minval(Exptrange_Mhmin_forSMdecays) to maxval(Exptrange_Mhmax_forSMdecays)
 ! otherwise this range is 
 !   Exptrange_Mhmin_forSMXS(x) to Exptrange_Mhmax_forSMXS(x)
 ! where x is the element of the array 'colliders'
  character(len=*),intent(in) :: str
  double precision,intent(in) :: mass
  logical :: inrange
  integer :: x, i
 
 !-TS(05/07/2012) Added a check whether the ranges are allocated (which happens in 
 !-subroutine setup_S95tables. If not, they are allocated and given default values.
 !-This is needed for HiggsSignals. 

  if(.not.allocated(Exptrange_Mhmin_forSMXS)) then
   allocate(Exptrange_Mhmin_forSMXS(size(colliders,dim=1)))
   allocate(Exptrange_Mhmax_forSMXS(size(colliders,dim=1)))
   do i=1, size(colliders)
    if(i.eq.get_collider_element_number('LEP ')) then
     Exptrange_Mhmin_forSMXS(i)=1.0D0
     Exptrange_Mhmax_forSMXS(i)=180.0D0
    else if(i.eq.get_collider_element_number('TEV ')) then
     Exptrange_Mhmin_forSMXS(i)=80.0D0
     Exptrange_Mhmax_forSMXS(i)=350.0D0
    else if(i.eq.get_collider_element_number('LHC7')) then
     Exptrange_Mhmin_forSMXS(i)=90.0D0
     Exptrange_Mhmax_forSMXS(i)=600.0D0
    else if(i.eq.get_collider_element_number('LHC8')) then
     Exptrange_Mhmin_forSMXS(i)=90.0D0
     Exptrange_Mhmax_forSMXS(i)=600.0D0
    else
     stop 'Error in subroutine inrange. str for XS unknown.'
    endif
   enddo 
  endif 
   
  if(.not.allocated(Exptrange_Mhmin_forSMdecays)) then
   allocate( Exptrange_Mhmin_forSMdecays(size(colliders,dim=1))) 
   allocate( Exptrange_Mhmax_forSMdecays(size(colliders,dim=1)))
   do i=1, size(colliders)
    if(i.eq.get_collider_element_number('LEP ')) then
     Exptrange_Mhmin_forSMdecays(i)=1.0D0
     Exptrange_Mhmax_forSMdecays(i)=180.0D0
    else if(i.eq.get_collider_element_number('TEV ')) then
     Exptrange_Mhmin_forSMdecays(i)=0.2D0
     Exptrange_Mhmax_forSMdecays(i)=350.0D0
    else if(i.eq.get_collider_element_number('LHC7')) then
     Exptrange_Mhmin_forSMdecays(i)=90.0D0
     Exptrange_Mhmax_forSMdecays(i)=600.0D0
    else if(i.eq.get_collider_element_number('LHC8')) then
     Exptrange_Mhmin_forSMdecays(i)=90.0D0
     Exptrange_Mhmax_forSMdecays(i)=600.0D0
    else
     stop 'Error in subroutine inrange. str for decay unknown.'
    endif
   enddo 
  endif
  
  if(str.eq.'SMBR')then   

    if(       (mass.gt.minval(Exptrange_Mhmin_forSMdecays)) &
         .and.(mass.lt.maxval(Exptrange_Mhmax_forSMdecays)) )then
      inrange=.True.
    else
      inrange=.False.
    endif 

  else
    x=get_collider_element_number(str)
 
    if(     (mass.gt.Exptrange_Mhmin_forSMXS(x)) &
       .and.(mass.lt.Exptrange_Mhmax_forSMXS(x)) )then
      inrange=.True.
    else
      inrange=.False.
    endif
   endif



 end function inrange   
 !**********************************************************
 function get_collider_element_number(collidername)
 ! this will return the position of the element with the value 'collidername' in the array 'colliders'
 !**********************************************************
  character(len=*),intent(in) :: collidername
  integer :: x,y
  integer :: get_collider_element_number

  y=0
  
!!    print *, collidername


  
  do x=lbound(colliders,dim=1),ubound(colliders,dim=1)
   if(index(colliders(x),collidername).gt.0)then
    get_collider_element_number=x
    y=y+1
   endif
  enddo

  if(y.ne.1)stop'problem in function get_collider_element_number'

 end function get_collider_element_number
 !****************************************************************** 
 function WhichColliderElement(expt,energy)
 ! this will return the position of the element corresponding to 'expt' in the array 'colliders'
 !****************************************************************** 
   use usefulbits, only : small
   character(LEN=3),intent(in) :: expt   
   integer :: WhichColliderElement
   double precision :: energy

   if(      expt.eq.'LEP' )then 
    WhichColliderElement=get_collider_element_number('LEP')
   elseif( (expt.eq.'CDF').or.(expt.eq.' D0').or.(expt.eq.'TCB') )then 
    WhichColliderElement=get_collider_element_number('TEV')
   elseif( (expt.eq.'ATL').or.(expt.eq.'CMS') )then 
    if(energy-7.0D0.le.small) then
     WhichColliderElement=get_collider_element_number('LHC7')
    else if(energy-8.0D0.le.small) then
     WhichColliderElement=get_collider_element_number('LHC8')
    else
     stop 'WhichColliderElement: Collider Energy not correctly specified.'
    endif 
   else
    WhichColliderElement=0
   endif

 end function WhichColliderElement
 !****************************************************************** 
 function WhichColliderString(expt, energy)
 ! this will return the contents of the element corresponding to 'expt' in the array 'colliders'
 !****************************************************************** 
   character(LEN=3),intent(in) :: expt   
   character(LEN=4) :: WhichColliderString
   double precision :: energy

   WhichColliderString=colliders(WhichColliderElement(expt, energy))

 end function WhichColliderString
 !******************************************************************      
 subroutine calcfact(proc,t,cfact,axis_i,axis_j,nc)
 !******************************************************************
 !for table type1, calls calcfact_t1
 !for table type2, calls calcfact_t2
 !****************************************************************** 
  use usefulbits, only : dataset,listprocesses  !internal
  implicit none      
  !--------------------------------------input
  type(dataset) :: t 
  type(listprocesses) :: proc
  !-----------------------------------output
  double precision :: cfact,axis_j,axis_i
  integer :: nc
  !-------------------------------------------      

  select case(proc%ttype)      
  case(1)
   call calcfact_t1(proc%tlist,proc%findj,t,cfact,axis_i,nc)
   axis_j=axis_i
  case(2)      
   call calcfact_t2(proc%tlist,proc%findj,proc%findi,t,cfact,axis_i,axis_j,nc)
  case default
! ### OS ### 
   cfact = -1D0
   return
   stop 'wrong input to function calcfact in module channels'
  end select            

 end subroutine calcfact

 !******************************************************************
 subroutine outputproc(proc,k,descrip,specific)
 !******************************************************************
 !for table type1, calls outputproc_t1
 !for table type2, calls outputproc_t2
 !if neither table applies, writes message
 !k is where output goes
 !specific=1 if the specific process should be printed 
 !e.g. ee->h1Z->bbZ
 !whereas specific==0 if the generic process should be printed
 !e.g. ee->hiZ->bbZ  
 !****************************************************************** 
  use usefulbits, only : listprocesses !input
  implicit none
  !--------------------------------------input
  integer,intent(in) :: k,specific
  integer :: i,j
  type(listprocesses),intent(in) :: proc
  character(LEN=200):: descrip  
  !-------------------------------------------         

  select case(specific)
  case(0)
   i=0
   j=0
  case(1)
   i=proc%findi
   j=proc%findj
  case default
  end select

  select case(proc%ttype)
  case(0)
   descrip='none of the processes apply'
  case(1)
   call outputproc_t1(proc%tlist,i,k,descrip)
  case(2)
   call outputproc_t2(proc%tlist,i,j,k,descrip) 
  case default
   
   stop 'wrong input to subroutine outputproc in module channels'
  end select              
            
 end subroutine outputproc 

 !******************************************************************      
 subroutine check_against_bound(proc,fact,axis_i,axis_j,ratio,predobs)
 !******************************************************************
 !for table type1, calls interpolate_tabletype1
 !for table type2, calls interpolate_tabletype2
 !****************************************************************** 
  use usefulbits, only : dataset,listprocesses  !internal
  use interpolate
  implicit none
  !--------------------------------------input
  integer :: predobs
  double precision :: fact,axis_i,axis_j
  type(listprocesses) :: proc  
  !-------------------------------------output
  double precision :: ratio  
  !-----------------------------------internal
  double precision :: Mi,Mj,interpol
  !-------------------------------------------      
  interpol=-1.0D0  
  Mi=axis_i
  Mj=axis_j            
  
  if(fact.gt.0.0D0)then 
   select case(proc%ttype)      
   case(1)
    call interpolate_tabletype1(Mi,S95_t1(proc%tlist),predobs,interpol)
   case(2)      
    call interpolate_tabletype2(Mi,Mj,S95_t2(proc%tlist),predobs,interpol)
   case default
! ### OS ###
     ratio = -1d0
      return
    write(*,*)'wrong input to subroutine check_against_bound in module channels'
    stop
   end select      
  endif      

!--TESTING!
!  if(predobs.eq.1) write(*,*) 'Interpolated value = ', interpol

  if(interpol.ge.0)then
   ratio=fact/interpol
  else
   ratio= -1.0D0
  endif

 end subroutine check_against_bound 

 !**********************************************************
 subroutine calcfact_t1(c,jj,t,cfact_t1,M_av,nc)
 !**********************************************************
 ! calculates fact for table type 1 
 ! Takes in to account how Standard Model-like parameter point is
 ! and whether there are any Higgs with slightly higher masses which
 ! can be combined with his result
 ! note: numerator and denominator are worked out separately
 !**********************************************************
  use usefulbits, only : dataset, np, div, vvsmall
  use theory_BRfunctions
  use theory_XS_SM_functions 
    
  implicit none      
  !--------------------------------------input
  type(dataset) :: t
  integer :: c,jj
  !-----------------------------------output
  double precision :: cfact_t1,M_av
  integer :: nc
  !-------------------------------------------
  integer :: f,j
  double precision :: M_tot
  double precision :: BR_Hbb_SM_av,BR_HWW_SM_av,BR_Htautau_SM_av
  double precision :: tev_XS_HW_SM_av,tev_XS_HZ_SM_av
  double precision :: tev_XS_H_SM_av
  double precision :: tev_XS_Hb_SM_av
  double precision :: tev_XS_ttH_SM_av
  double precision :: lhc7_XS_H_SM_av
  double precision :: lhc7_XS_VBF_SM_av
  double precision :: BR_Zll,BR_Znunu,BR_Wlnu,BR_Ztautau
  double precision :: BR_Whad,BR_Zhad  
  double precision,allocatable :: mass(:),fact(:)
  integer,allocatable :: model_like(:)
  integer :: npart !number of particles

  !source: PDG (Yao et al. J Phys G 33 (2006)) 
  BR_Zll=3.363D-2+3.366D-2      !BR_Zll = sum(l=e,mu), BR(Z ->l+ l-)
  BR_Znunu=20D-2                !BR_Znunu = BR(Z ->nu_l nu_l-bar) ('invisible')
  BR_Wlnu=10.75D-2+10.57D-2     !BR_Wlnu = sum(l=e,mu),
                                !BR(W+ ->l+ nu_l) = BR(W- ->l- + nu_l-bar)
  BR_Whad=67.6D-2     
  BR_Zhad=69.91D-2         

  BR_Ztautau=3.370D-2

  npart=np( S95_t1(c)%particle_x )
  allocate(mass(npart),fact(npart),model_like(npart))

  mass(:)=t%particle( S95_t1(c)%particle_x )%M(:) 

  fact= 0.0D0     
  cfact_t1=0.0D0    
  model_like=0

  !now calculate numerator of 'fact' 
  do j=1,npart    
  
   if(     (abs(mass(jj)-mass(j)).le.S95_t1(c)%deltax) &
    & .and.(            mass(jj).le.mass(j)           )  )then
         
    select case(S95_t1(c)%id) !these can be compacted, but not doing it yet... doing it at the 
                              !same time as revamping the model_likeness test
    case(142)
     fact(j)=t%lep%XS_hjZ_ratio(j)   *t%BR_hjbb(j) !notice that this is not absolute XS
    case(143)
     fact(j)=t%lep%XS_hjZ_ratio(j)   *t%BR_hjtautau(j) !notice that this is not absolute XS
    case(300)
     fact(j)=t%lep%XS_hjZ_ratio(j)   !notice that this is not absolute XS
    case(400,401,402,403)
     fact(j)=t%lep%XS_hjZ_ratio(j)   *t%BR_hjinvisible(j)   !notice that this is not absolute XS
    case(500)
     fact(j)=t%lep%XS_hjZ_ratio(j)   *t%BR_hjgaga(j)   !notice that this is not absolute XS
    case(600)
     fact(j)=t%lep%XS_hjZ_ratio(j)   &!notice that this is not absolute XS
      & *(t%BR_hjss(j)+t%BR_hjcc(j)+t%BR_hjbb(j)+t%BR_hjgg(j))   
    case(711,713)
     call model_likeness(j,S95_t1(c)%id,t,model_like(j),fact(j))
    case(721,723,741,743)
     call model_likeness(j,S95_t1(c)%id,t,model_like(j),fact(j))
    case(731,733)
     call model_likeness(j,S95_t1(c)%id,t,model_like(j),fact(j))
    case(801,811,821)
     fact(j)=t%lep%XS_HpjHmj_ratio(j)*(t%BR_Hpjcs(j)+t%BR_Hpjcb(j))**2.0D0 !notice that this is not absolute XS      
    case(802)
     fact(j)=t%lep%XS_HpjHmj_ratio(j)*(t%BR_Hpjcs(j)+t%BR_Hpjcb(j))*2.0D0*t%BR_Hpjtaunu(j) !notice that this is not absolute XS       
    case(803,813)
     fact(j)=t%lep%XS_HpjHmj_ratio(j)*t%BR_Hpjtaunu(j)**2.0D0 !notice that this is not absolute XS     
    case(8742,4493,9475,5482,5570,5876,1024,9889,3534,6089,10235,3047,10799,3564,6166,6296)   
     fact(j)=t%tev%XS_hjZ_ratio(j)   *t%tev%XS_HZ_SM(j)  *t%BR_hjbb(j)        
    case(8958,5489,5624,9236,3930,6039,3216,10433,6221,10600)
     fact(j)=t%tev%XS_hj_ratio(j)*t%tev%XS_H_SM(j)   *t%BR_hjWW(j)
    case(3331)
     call model_likeness(j,S95_t1(c)%id,t,model_like(j),fact(j))
     fact(j)=fact(j)*t%tev%XS_H_SM(j)*t%BR_HWW_SM(j)!to get the normalisation right
    case(5757)
     call model_likeness(j,S95_t1(c)%id,t,model_like(j),fact(j))
    case(8957,5472,9219,9463,9596,5828,1970,3493,5972,3155,5613,9868,10068,6092,10217,10239,10796,0874,6220)  
     fact(j)=t%tev%XS_hjW_ratio(j)   *t%tev%XS_HW_SM(j)  *t%BR_hjbb(j)  
    case(5485,7307,5873)   
     fact(j)=t%tev%XS_hjW_ratio(j)   *t%tev%XS_HW_SM(j)  *t%BR_hjWW(j)       
    case(9071,2491,5740,5980,1014,3363,4555)
     fact(j)=t%tev%XS_hj_ratio(j)    *t%tev%XS_H_SM(j)   *t%BR_hjtautau(j) 
    case(8961,0598)   
     call model_likeness(j,S95_t1(c)%id,t,model_like(j),fact(j))    
    case(9284,5503,5726,10105)  
     fact(j)=t%tev%XS_hjb_ratio(j)    *t%tev%XS_Hb_SM(j)  *t%BR_hjbb(j)/0.9D0/2.0D0 
    case(6083)  
     fact(j)=t%tev%XS_hjb_ratio(j)    *t%tev%XS_Hb_SM(j)  *t%BR_hjtautau(j)/0.1D0/2.0D0 
    case(1514,5601,5737)  
     fact(j)=(  t%tev%XS_hjZ_ratio(j) * t%tev%XS_HZ_SM(j)      &
      &  +   t%tev%XS_hjW_ratio(j) * t%tev%XS_HW_SM(j)      &
      &  +   t%tev%XS_hj_ratio(j)  * t%tev%XS_H_SM(j)       &
      &  +   t%tev%XS_vbf_ratio(j) * t%tev%XS_vbf_SM(j)     &
      &    )                  & 
      &  *   t%BR_hjgaga(j)                    
    case(5858,6177,6295,1887,10065,10485,4960)  
     call model_likeness(j,S95_t1(c)%id,t,model_like(j),fact(j)) 
    case(7081,9166,9483,5586,9642,1266,0432,9891,5285,3935,6087,6170,10212,6223,6299,10583,10798,10596)
     call model_likeness(j,S95_t1(c)%id,t,model_like(j),fact(j))
    case(2012015)
     call model_likeness(j,S95_t1(c)%id,t,model_like(j),fact(j))
    case(10010,10607,6436)
     call model_likeness(j,S95_t1(c)%id,t,model_like(j),fact(j))  
    case(9248,10133,10439)
     call model_likeness(j,S95_t1(c)%id,t,model_like(j),fact(j))    
    case(4800)
     call model_likeness(j,S95_t1(c)%id,t,model_like(j),fact(j))    
    case(5845)
     call model_likeness(j,S95_t1(c)%id,t,model_like(j),fact(j))
    case(6171)
     call model_likeness(j,S95_t1(c)%id,t,model_like(j),fact(j))     
    case(9290)
     call model_likeness(j,S95_t1(c)%id,t,model_like(j),fact(j)) 
    case(5984,9714,6006,4481,6095,10432,6179,6302,10599)
     call model_likeness(j,S95_t1(c)%id,t,model_like(j),fact(j))
    case(3556,1931)
     fact(j)=t%tev%XS_hjb_ratio(j)     *t%tev%XS_Hb_c2_SM(j)   *t%BR_hjbb(j) 
    case(4782)
     fact(j)=t%tev%XS_hjb_ratio(j)     *t%tev%XS_Hb_c1_SM(j)   *t%BR_hjbb(j)                                        
    case(9465,5871,9022,0710,9887,4162,10102,4468)
     call model_likeness(j,S95_t1(c)%id,t,model_like(j),fact(j))
    !case(9023)
    ! call model_likeness(j,S95_t1(c)%id,t,model_like(j),fact(j))
    case(9674)
     call model_likeness(j,S95_t1(c)%id,t,model_like(j),fact(j))
    case(9897,9999)
     call model_likeness(j,S95_t1(c)%id,t,model_like(j),fact(j))
    case(0024,5985,0968,5974,4885)
     fact(j)=t%tev%XS_hjb_ratio(j)     *t%tev%XS_Hb_c3_SM(j)   *t%BR_hjtautau(j) 
    case(5739,10574)
     call model_likeness(j,S95_t1(c)%id,t,model_like(j),fact(j))
     fact(j)=fact(j)*t%tev%XS_ttH_SM(j)*t%BR_Hbb_SM(j)!to get the normalisation right
    case(2012135,12025)
     call model_likeness(j,S95_t1(c)%id,t,model_like(j),fact(j))
    case(0611)
     fact(j)=t%tev%XS_hj_ratio(j)*t%tev%XS_H_SM(j)   *t%BR_hjZga(j)
    case(1269,1270,2011094)
     call model_likeness(j,S95_t1(c)%id,t,model_like(j),fact(j))
    case(1811)
     call model_likeness(j,S95_t1(c)%id,t,model_like(j),fact(j))
    case(1812,2011138,2011151,2760,2013090,8353,7712,11002,11008)
     call model_likeness(j,S95_t1(c)%id,t,model_like(j),fact(j))
    case(6008,9998)
     call model_likeness(j,S95_t1(c)%id,t,model_like(j),fact(j))
    case(6082,6182,6219,6276,10573)
     call model_likeness(j,S95_t1(c)%id,t,model_like(j),fact(j))
    case(6096,10606,10806,10884)
     call model_likeness(j,S95_t1(c)%id,t,model_like(j),fact(j))
    case(6183,3233)
     call model_likeness(j,S95_t1(c)%id,t,model_like(j),fact(j))  
    case(6229)
     call model_likeness(j,S95_t1(c)%id,t,model_like(j),fact(j))  
    case(6286,6301,6304,6305,6309)
     call model_likeness(j,S95_t1(c)%id,t,model_like(j),fact(j)) 
    case(6091,1268)
     call model_likeness(j,S95_t1(c)%id,t,model_like(j),fact(j))
    case(2011048,11004,11015,11013,11028,11006,11017,110271,110272,14161,14162,5064,2012017,2011150,2011131)
     call model_likeness(j,S95_t1(c)%id,t,model_like(j),fact(j))
    case(2011162,1415,2012092,20130131)
     call model_likeness(j,S95_t1(c)%id,t,model_like(j),fact(j))
    case(20130132)
     fact(j)=t%lhc8%XS_hj_ratio(j)*t%lhc8%XS_H_SM(j)*t%BR_hjZZ(j)
    case(20130133)
     fact(j)=(t%lhc8%XS_hjZ_ratio(j) * t%lhc8%XS_HZ_SM(j)     &
      &  +   t%lhc8%XS_hjW_ratio(j) * t%lhc8%XS_HW_SM(j)      &
      &  +   t%lhc8%XS_vbf_ratio(j) * t%lhc8%XS_vbf_SM(j))    & 
      &  *   t%BR_hjZZ(j)     
    case(11025,1997,12041,130021,130022)
     call model_likeness(j,S95_t1(c)%id,t,model_like(j),fact(j))
    case(2011026,11005,11016,11026,3478,3357,2011148,2012016)
     call model_likeness(j,S95_t1(c)%id,t,model_like(j),fact(j))
    case(110212)
     call model_likeness(j,S95_t1(c)%id,t,model_like(j),fact(j))
    case(2011025,2011085,2011161,5895,1414,2012091,2012168,1487,12001,12015,13001,11010,11030,11021)
     call model_likeness(j,S95_t1(c)%id,t,model_like(j),fact(j))
    case(13006, 13075515,2013009)
     call model_likeness(j,S95_t1(c)%id,t,model_like(j),fact(j))
    case(11031,12044,13012)
     call model_likeness(j,S95_t1(c)%id,t,model_like(j),fact(j))     
    case(13011)
     call model_likeness(j,S95_t1(c)%id,t,model_like(j),fact(j))     
    case(11034,12039,13009,2012078,12006,12051)
     call model_likeness(j,S95_t1(c)%id,t,model_like(j),fact(j))     
    case(2011005,3615,2012018,12046)
     call model_likeness(j,S95_t1(c)%id,t,model_like(j),fact(j))
    case(2748, 1408, 2012019)
     call model_likeness(j,S95_t1(c)%id,t,model_like(j),fact(j))
    case(7214)
     call model_likeness(j,S95_t1(c)%id,t,model_like(j),fact(j))
    case(10500)
     call model_likeness(j,S95_t1(c)%id,t,model_like(j),fact(j))
    case(11003,11014,2577,11024,1489,12042,13003)
     call model_likeness(j,S95_t1(c)%id,t,model_like(j),fact(j))
    case(2011020,2011021)
     fact(j)=t%lhc7%XS_hj_ratio(j)*t%lhc7%XS_H_SM(j)   *t%BR_hjmumu(j)
    case(5429,2011052,2011111,2011134)
     fact(j)=t%lhc7%XS_hj_ratio(j)*t%lhc7%XS_H_SM(j)   *t%BR_hjWW(j)
    case(2012012,2012158,2013030)
     call model_likeness(j,S95_t1(c)%id,t,model_like(j),fact(j))
    case(10002,5003,2011132,2012094,110201,110292)
     fact(j)=t%lhc7%XS_hj_ratio(j)*t%lhc7%XS_H_SM(j)   *t%BR_hjtautau(j)
    case(12050)
     fact(j)=t%lhc8%XS_hj_ratio(j)*t%lhc8%XS_H_SM(j)   *t%BR_hjtautau(j)
    case(11009,11020,2011133,2012014,2012160)
     call model_likeness(j,S95_t1(c)%id,t,model_like(j),fact(j))
    case(110291,12043)
     call model_likeness(j,S95_t1(c)%id,t,model_like(j),fact(j))
    case(2013010)
     call model_likeness(j,S95_t1(c)%id,t,model_like(j),fact(j))
    case(11011,2011163,11022,11032,1488,12008,12045,2011157)
     call model_likeness(j,S95_t1(c)%id,t,model_like(j),fact(j))
    case(2011103,2012161,11012)
     call model_likeness(j,S95_t1(c)%id,t,model_like(j),fact(j))
    case(13022)
     fact(j)=t%lhc8%XS_vbf_ratio(j)*div(t%BR_hjWW(j),t%BR_HWW_SM(j) ,0.0D0,1.0D0)
    case(13013)
     fact(j)=t%lhc8%XS_vbf_ratio(j)*t%BR_hjinvisible(j)
    case(13018)
     fact(j)=t%lhc8%XS_hjZ_ratio(j)*t%BR_hjinvisible(j)
    case(2013011)
!    Limit is on sigma(HZ)*BR(H->inv)*(BR(Z->ll)+BR(Z->tautau) 
!    Data given in fb - (multiply by 1000)
          fact(j)=1000.D0*t%lhc8%XS_hjZ_ratio(j)*t%lhc8%XS_HZ_SM(j)   &    
      &    *t%BR_hjinvisible(j)*(BR_Zll+BR_Ztautau)
            
!       print *, 1000.D0*t%lhc8%XS_hjZ_ratio(j)*t%lhc8%XS_HZ_SM(j)*(BR_Zll+BR_Ztautau), fact(j)
      
    case(2011112)
     call model_likeness(j,S95_t1(c)%id,t,model_like(j),fact(j))
    case(2011135)
     call model_likeness(j,S95_t1(c)%id,t,model_like(j),fact(j))
    case(6224,6225,6226)
     call model_likeness(j,S95_t1(c)%id,t,model_like(j),fact(j))
    case default
     stop 'wrong input to function calcfact_t1 in module S95tables'
    end select

   endif

  enddo
 
  if(fact(jj).le.vvsmall)then!A  !Higgs jj doesn't contribute - wait until another call of this subroutine before
                               !looking at nearby masses
   M_av = mass(jj)
   nc=0
   cfact_t1=0.0D0   
  else!A
   !find M_av (only using higgs which have non-zero fact): 
   f=0  
   M_tot=0.0D0
   do j=1,npart    
    if(   fact(j).gt.vvsmall  )then       
     f=f+1
     M_tot=M_tot+mass(j) 
    endif
   enddo 
   
   nc=f !f will always be > 0 because we've already made sure that fact(jj)>0.0D0
  
   M_av = M_tot/dble(nc)  
   
   if((WhichColliderString(S95_t1(c)%expt,S95_t1(c)%energy).eq.'LEP'))then!B 
    cfact_t1=sum(fact)
   elseif(S95_t1(c)%particle_x .ne. Hneut)then!B
    cfact_t1=sum(fact)
   else!B

    if(f.eq.1)then !have already calculated these in theo_manip to save time
     BR_Hbb_SM_av     = t%BR_Hbb_SM(jj)
     BR_HWW_SM_av     = t%BR_HWW_SM(jj)
     BR_Htautau_SM_av = t%BR_Htautau_SM(jj)
     tev_XS_HW_SM_av  = t%tev%XS_HW_SM(jj)
     tev_XS_HZ_SM_av  = t%tev%XS_HZ_SM(jj)
     tev_XS_H_SM_av   = t%tev%XS_H_SM(jj)
     tev_XS_Hb_SM_av  = t%tev%XS_Hb_SM(jj)
     tev_XS_ttH_SM_av = t%tev%XS_ttH_SM(jj)
     lhc7_XS_H_SM_av  = t%lhc7%XS_H_SM(jj)
     lhc7_XS_VBF_SM_av= t%lhc7%XS_vbf_SM(jj)
    else   
     BR_Hbb_SM_av     = BRSM_Hbb(M_av)  
     BR_HWW_SM_av     = BRSM_HWW(M_av)
     BR_Htautau_SM_av = BRSM_Htautau(M_av)
     tev_XS_HW_SM_av  = XS_tev_HW_SM(M_av)
     tev_XS_HZ_SM_av  = XS_tev_HZ_SM(M_av)
     tev_XS_H_SM_av   = XS_tev_gg_H_SM(M_av)+XS_tev_bb_H_SM(M_av) 
     tev_XS_Hb_SM_av  = XS_tev_bg_Hb_SM(M_av) 
     tev_XS_ttH_SM_av = XS_tev_ttH_SM(M_av)    
     lhc7_XS_H_SM_av  = XS_lhc7_gg_H_SM(M_av)+XS_lhc7_bb_H_SM(M_av) 
     lhc7_XS_VBF_SM_av= XS_lhc7_vbf_SM(M_av)
    endif     
    
    ! now include denominator of 'fact'
    select case(S95_t1(c)%id)
    case(8742,5482,5570,4493,9475,5876,1024,9889,3534,6089,10235,3047,10799,3564,6166,6296)   
     do j=1,npart
      fact(j)=    div( fact(j) ,  tev_XS_HZ_SM_av  *  BR_Hbb_SM_av    ,0.0D0,0.0D0) 
     enddo            
    case(8958,5489,5624,9236,3930)
     do j=1,npart
      fact(j)= div( fact(j) ,  tev_XS_H_SM_av   *  BR_HWW_SM_av    ,0.0D0,0.0D0)  
     enddo     
    case(8957,5472,9219,9463,9596,5828,1970,3493,5972,3155,5613,9868,10068,6092,10217,10239,10796,0874,6220) 
     do j=1,npart 
      fact(j)=    div( fact(j) ,  tev_XS_HW_SM_av  *  BR_Hbb_SM_av    ,0.0D0,0.0D0) 
     enddo           
    case(5503,9284,5726,10105)  
     do j=1,npart
      fact(j)=    div( fact(j) ,  tev_XS_Hb_SM_av  ,0.0D0,0.0D0) 
     enddo            
    case(6083)  
     do j=1,npart
      fact(j)=    div( fact(j) ,  tev_XS_Hb_SM_av ,0.0D0,0.0D0) 
     enddo           
    case(7307,5873) 
     do j=1,npart
      fact(j)= div( fact(j) ,  tev_XS_HW_SM_av  *  BR_HWW_SM_av    ,0.0D0,0.0D0)                                             
     enddo                         
    case(8961,0598,10010,9290,9674,9897,9999,10607,6436)
    case(7081,9166,9483,5586,9642,1266,0432,9891,5285,3935,6087,6170,10212,6223,6299,10583,10798,10596)
    case(2012015)
    case(9248,5845,4800,5858,6177,6295,1887,10065,10485,10133,10439,4960) 
    case(9465,5871,9022,9023,0710,9887)
    case(6171)
    case(6183,3233) 
    case(6229)    
    case(6304)    
    case(5984,9714,4162,10102,6006,4481,4468,5757,6095,10432,6179,6302,10599)
    case(5739,10574)
     do j=1,npart
      fact(j)=    div( fact(j) ,  tev_XS_ttH_SM_av * BR_Hbb_SM_av     ,0.0D0,0.0D0) 
     enddo         
    case(2012135,12025)   
    case(6008,9998)
    case(6082,6182,6219,6276,10573)
    case(6096,10606,10806,10884)
    case(6091,1268)
    case(5485,9071,2491,5601,1514,3556,5740,4555,0024,5980,1014,5985, &
      &  0611,3363,6039,3216,0968,5974,1931,4885,6221)
    case(3331)
    case(6286,6301,6305,6309)
    case(10600)
    case(10433)
     do j=1,npart
      fact(j)=div(fact(j)  ,  tev_XS_H_SM_av   *  BR_HWW_SM_av     ,0.0D0,0.0D0)  
     enddo  
    case(2011048,11004,11015,11013,11028,11006,11017,110271,110272,14161,14162,5064,2012017,2011150,2011131)
    case(2011162,1415,2012092,20130131)
    case(11025,1997,12041,130021,130022)
    case(2011026,11005,11016,11026,3478,3357,2011148,2012016)
    case(110212)
    case(2011025,2011085,2011161,5895,1414,2012091,2012168,1487,12001,12015,13001,11010,11030,11021)
    case(13006, 13075515,2013009)
    case(11031,12044,13012)   
    case(13011)   
    case(11034,12039,13009,2012078,12006,12051)    
    case(2011005,3615,2012018,12046) 
    case(2748, 1408, 2012019) 
    case(7214)
    case(4782)
    case(5429,2011052,2011111,2011134) 
     do j=1,npart
      fact(j)= div( fact(j) ,  lhc7_XS_H_SM_av   *  BR_HWW_SM_av    ,0.0D0,0.0D0)  
     enddo     
    case(2012012,2012158,2013030)
    case(10002,5003,2011132,2012094,110201,110292,12050) 
    case(20130132,20130133)
    case(11009,11020,2011133,2012014,2012160)
    case(110291,12043)
    case(2013010)
    case(11002,11008) 
    case(11003,11014,2577,11024,1489,12042,13003)
    case(2011020,2011021)
    case(10500)   
    case(11011,2011163,11022,11032,1488,12008,12045,2011157)   
    case(2011103,2012161,11012) 
    case(13022)
    case(13013)
    case(13018)
    case(2013011)
    case(2011112)
    case(2011135) 
    case(6224,6225,6226)
    case default
     stop 'error calculating denom. in calcfact_t1'
    end select
   
    cfact_t1=sum(fact)

   endif!B
  endif!A

  deallocate(mass)
  deallocate(fact)
  deallocate(model_like)

 end subroutine calcfact_t1  
 !**********************************************************      
 subroutine calcfact_t2(c,jj,ii,t,cfact_t2,axis_i,axis_j,nc)
 !**********************************************************
 !calculates fact for table type 2 
 !**********************************************************      
  use usefulbits, only : dataset,np,vsmall,not_a_particle
  implicit none      
  !--------------------------------------input      
  type(dataset) :: t
  integer :: c,jj,ii
  !-----------------------------------output
  double precision :: cfact_t2,axis_i,axis_j
  integer :: nc
  !-------------------------------------------
  integer :: f,i,j,npart2,npart1
  double precision :: fact,eps2,crosssection,Mi_av,Mj_av,masstot
  double precision,allocatable :: massj(:),massi(:)

   eps2=0.02D0

   npart2=np( S95_t2(c)%particle_x2 )
   allocate(massj(npart2))
   massj(:)=t%particle( S95_t2(c)%particle_x2 )%M(:)

   if(S95_t2(c)%particle_x1.eq.not_a_particle)then
     npart1=1
     allocate(massi(npart1))
     massi(:)=-1.0D0
   else
     npart1=np( S95_t2(c)%particle_x1 )
     allocate(massi(npart1))
     massi(:)=t%particle( S95_t2(c)%particle_x1 )%M(:)
   endif

   Mj_av=massj(jj)
   Mi_av=massi(ii)

   fact= 0.0D0                    
   cfact_t2=0.0D0     
   masstot=0.0D0     
   j=jj
   i=ii
   f=1
        
   select case(S95_t2(c)%id)
   case(150)     
    fact=test_appl(t%lep%XS_hjZ_ratio(j)*t%BR_hjhihi(j,i)*t%BR_hjbb(i)**2.0D0)!table 15 hep-ex/0602042 XS ratio
   case(160)  
    fact=test_appl(t%lep%XS_hjZ_ratio(j)*t%BR_hjhihi(j,i)*t%BR_hjtautau(i)**2.0D0)!table 16 hep-ex/0602042 XS ratio
   case(180)  
    fact=test_appl(t%lep%XS_hjhi_ratio(j,i)*t%BR_hjbb(j)*t%BR_hjbb(i))!table 18 hep-ex/0602042 XS ratio         
   case(190)
    fact=test_appl(t%lep%XS_hjhi_ratio(j,i)*t%BR_hjtautau(j)*t%BR_hjtautau(i))!table 19 hep-ex/0602042 XS ratio
   case(200)  
    fact=test_appl(t%lep%XS_hjhi_ratio(j,i)*t%BR_hjhihi(j,i)*t%BR_hjbb(i)**3.0D0)!table 20 hep-ex/0602042 XS ratio
   case(210) 
    fact=test_appl(t%lep%XS_hjhi_ratio(j,i)*t%BR_hjhihi(j,i)*t%BR_hjtautau(i)**3.0D0)!table 21 hep-ex/0602042 XS ratio
   case(220) 
    fact=test_appl(t%lep%XS_hjZ_ratio(j)*t%BR_hjhihi(j,i)*t%BR_hjbb(i)*t%BR_hjtautau(i))!table 22 hep-ex/0602042 XS ratio
   case(230)                     
    fact=test_appl(t%lep%XS_hjhi_ratio(j,i)*t%BR_hjbb(j)*t%BR_hjtautau(i))!table 23 hep-ex/0602042 XS ratio        
   case(240) 
    fact=test_appl(t%lep%XS_hjhi_ratio(j,i)*t%BR_hjtautau(j)*t%BR_hjbb(i))!table 24 hep-ex/0602042 XS ratio
   case(905) 
    fact=test_appl(t%lep%XS_CpjCmj(j)*t%BR_CjqqNi(j,i)**2.0D0)!fig 5 hep-ex/0401026 absolute XS in fb
   case(906) 
    fact=test_appl(t%lep%XS_CpjCmj(j)*t%BR_CjqqNi(j,i)*t%BR_CjlnuNi(j,i))!fig 6 hep-ex/0401026 absolute XS in fb 
   case(907) 
    fact=test_appl( t%lep%XS_CpjCmj(j)*t%BR_CjlnuNi(j,i)**2.0D0)!fig 7 hep-ex/0401026 absolute XS in fb
   case(908) 
    fact=test_appl(t%lep%XS_CpjCmj(j)) !fig 8 hep-ex/0401026 absolute XS in fb
   case(909)
    fact=test_appl( t%lep%XS_NjNi(j,i)*t%BR_NjqqNi(j,i)) !fig 9 hep-ex/0401026 absolute XS in fb
   case(910)
    fact=test_appl(t%lep%XS_NjNi(j,i))!fig 10 hep-ex/0401026 absolute XS in fb
   case(3381)
    fact=test_appl(t%tev%XS_hj_ratio(j) * t%tev%XS_H_SM(j) * t%BR_hjhihi(j,i) * &
                & t%BR_hjmumu(i)**2.0D0 )! arXiv:0905.3381 table I, absolute XS in fb
   case(3382)
    fact=test_appl(t%tev%XS_hj_ratio(j) * t%tev%XS_H_SM(j) * t%BR_hjhihi(j,i) * &
                 & 2.0D0 * t%BR_hjtautau(i) * t%BR_hjmumu(i) )! arXiv:0905.3381 table II (also using fig 3b), absolute XS in fb
   case(6227)

    f=0
    do j=1,npart2
     if(     (abs(massj(jj)-massj(j)).le.S95_t2(c)%deltax) &
       & .and.(            massj(jj).le.massj(j)           )  )then

      crosssection=test_appl( t%tev%XS_hjb_ratio(j)*t%tev%XS_Hb_c4_SM(j) )

      if(crosssection.gt.vsmall)then
        f=f+1
        fact=fact+crosssection
        masstot=massj(j)+masstot
      endif

     endif
    enddo
    
    if(f.ne.0)then
      Mj_av=masstot/dble(f)
    endif

   case default
    stop 'wrong input to function calcfact_t2 in module S95tables'
   end select  

   if(S95_t2(c)%particle_x1.eq.not_a_particle)then
     select case(S95_t2(c)%id)
     case(6227)
     axis_i=t%BR_hjtautau(jj) 
     case default 
      stop'Problem in subroutine calcfact_t2 (y1)'
     end select
   else
     axis_i=Mi_av
   endif

   if(S95_t2(c)%particle_x2.eq.not_a_particle)then
     select case(S95_t2(c)%id)
     case default 
       stop'Problem in subroutine calcfact_t2 (y2)'
     end select
   else
     axis_j=Mj_av
   endif

   cfact_t2=cfact_t2+fact      
  
   nc=f

   deallocate(massi)
   deallocate(massj)
 
  contains
  
  !********************************************************      
  function test_appl(x)  
  !********************************************************
   implicit none
   !--------------------------------------input
   double precision :: x
   !-----------------------------------function
   double precision :: test_appl
   !-------------------------------------------   
         
    select case(S95_t2(c)%id)
    case(150,160,180,190,200,210,220,230,240,3381,3382)
     if(S95_t2(c)%needs_M2_gt_2M1.and.(massj(j).lt.2.0D0*massi(i)))then
      test_appl=0.0D0 
     elseif(massj(j).lt.massi(i))then
      test_appl=0.0D0      
     else
      test_appl=x
     endif
    case(905,906,907,909)
     if(abs(minval(massi)-massi(i)).gt.vsmall)then !checking that lightest neutralino in process is lightest neutralino in model
      test_appl=0.0D0
     elseif(massj(j).lt.massi(i))then
      test_appl=0.0D0   
     else 
      test_appl=x 
     endif
    case(908)
     if( abs(t%BR_CjWNi(j,i)-1.0D0) .gt. eps2 )then
      test_appl=0.0D0
     elseif(abs(minval(massi)-massi(i)).gt.vsmall)then !checking that lightest neutralino in process is lightest neutralino in model
      test_appl=0.0D0
     elseif(massj(j).lt.massi(i))then
      test_appl=0.0D0   
     else
      test_appl=x 
     endif
    case(910)
     if( abs(t%BR_NjZNi(j,i)-1.0D0) .gt. eps2 )then
      test_appl=0.0D0
     elseif(abs(minval(massi)-massi(i)).gt.vsmall)then !checking that lightest neutralino in process is lightest neutralino in model
      test_appl=0.0D0
     elseif(massj(j).lt.massi(i))then
      test_appl=0.0D0   
     else
      test_appl=x 
     endif
    case(6227)
     if( ( t%BR_hjtautau(j)+t%BR_hjbb(j) ).le.0.98D0)then
      test_appl=0.0D0   
     else
      test_appl=x 
     endif
    case default
      stop'error in function test_appl'
    end select              

  end function test_appl
      
 end subroutine calcfact_t2            
      
 !********************************************************      
 subroutine outputproc_t1(tlistn,jj,k,descrip)
 !********************************************************            
 ! uses information about the process to output a description
 ! for processes using table type 1
 ! note: at the moment, np(x) (and so ii and jj) needs to be 1 digit long i.e. nH<10
 !******************************************************** 
  implicit none
  !--------------------------------------input
  integer :: tlistn
  integer :: jj,k  
  !-----------------------------------internal
  character(LEN=1) :: j
  character(LEN=45) :: label
  character(LEN=200):: descrip
  !-------------------------------------------
         
  if(jj.ne.0)then
   write(j,'(I1)')jj      
  else
   j='j'
  endif

  if(k.eq.21)then
   label=''  !no need to lable each line in Key.dat
  else
   label='('//trim(S95_t1(tlistn)%label)//')'
  endif          
  
  descrip=''   
  
  select case(S95_t1(tlistn)%id)
  case(142)
   descrip=' (e e)->(h'//j//')Z->(b b-bar)Z   '     //label             
  case(143)
   descrip=' (e e)->(h'//j//')Z->(tau tau)Z   ' //label
  case(300)
   descrip=' (e e)->(h'//j//')Z->(...)Z   '     //label 
  case(400,401,402,403)
   descrip=' (e e)->(h'//j//')Z->(invisible)Z   '     //label   
  case(500)
   descrip=' (e e)->(h'//j//')Z->(gamma gamma)Z   '     //label  
  case(600)
   descrip=' (e e)->(h'//j//')Z->(2 jets)Z   '     //label 
  case(711)
   descrip=' (e e)->b b-bar(h'//j//')->b b-bar(b b-bar) where h'//j//' is CP even '     //label 
  case(713)
   descrip=' (e e)->b b-bar(h'//j//')->b b-bar(b b-bar) where h'//j//' is CP odd '     //label 
  case(721,741)
   descrip=' (e e)->b b-bar(h'//j//')->b b-bar(tau tau) where h'//j//' is CP even '     //label 
  case(723,743)
   descrip=' (e e)->b b-bar(h'//j//')->b b-bar(tau tau) where h'//j//' is CP odd '     //label  
  case(731)
   descrip=' (e e)->tau tau(h'//j//')->tau tau(tau tau) where h'//j//' is CP even '     //label 
  case(733)
   descrip=' (e e)->tau tau(h'//j//')->tau tau(tau tau) where h'//j//' is CP odd '     //label  
  case(801,811,821)
   descrip=' (e e)->(H'//j//'+)(H'//j//'-)->4 quarks  '     //label
  case(802)
   descrip=' (e e)->(H'//j//'+)(H'//j//'-)->(2 quarks) tau nu '     //label
  case(803,813)
   descrip=' (e e)->(H'//j//'+)(H'//j//'-)->tau nu tau nu'     //label
  case(5482,5570,8742,4493,9475,5876,1024,9889,3534,6089,10235,3047,10799,3564,6166,6296)  
   descrip=' (p p-bar)->Z(h'//j//')->l l (b b-bar)   '  //label  
  case(9236,3930,8958,6039,3216,10433,6221,10600)  
   descrip=' (p p-bar)->h'//j//'->W W   '           //label    
  case(9219,9463,5472,8957,9596,5828,1970,3493,5972,3155,5613,9868,10068,6092,10217,10239,10796,0874,6220,6309)  
   descrip=' (p p-bar)->W(h'//j//')->l nu (b b-bar)   ' //label     
  case(5489)  
   descrip=' (p p-bar)->h'//j//'->W W->e mu   '     //label     
  case(5624)  
   descrip=' (p p-bar)->h'//j//'->W W->l l   '      //label
  case(3331)  
   descrip=' (p p-bar)->h'//j//'->V V   '      //label
  case(5757)  
   descrip=' (p p-bar)->h'//j//'/VBF->W W->l l where h'//j//' is SM-like   '      //label
  case(5485,5873)  
   descrip=' (p p-bar)->W(h'//j//')->W W W->l l nu nu   '         //label  
  case(9071,2491,5740,5980,1014,3363,4555)  
   descrip=' (p p-bar)->h'//j//'->tau tau   '       //label         
  case(8961,9465,9290,9713,9674,0598,9897,9998,9999,6008,6096,6183,3233,6229,6304,10606,10806,10884)  
   descrip=' (p p-bar)->h'//j//'+... where h'//j//' is SM-like  ' //label        
  case(9284,5503,5726,3556,10105,1931,4782)  
   descrip=' (p p-bar)->h'//j//'(b/b-bar)->(b b-bar) (b/b-bar)   '          //label  
  case(6224,6225,6226)  
   descrip=' (p p-bar)->h'//j//'(b/b-bar)->(b b-bar) (b/b-bar) or (tau tau) (b/b-bar) '     //label 
  case(7307)  
   descrip=' (p p-bar)->W(h'//j//')->W W W   '      //label 
  case(6301)  
   descrip=' (p p-bar)->V h'//j//'->V W W   '      //label 
  case(5601,5737,1514)    
   descrip=' (p p-bar)->h'//j//'+...->gamma gamma+... '  //label 
  case(5858,6177,6295,1887,10065,10485,4960)    
   descrip=' (p p-bar)->h'//j//'+...->gamma gamma+... where h'//j//' is SM-like  '  //label 
  case(7081,9166,9483,5586,9642,1266,0432,9891,5285,3935,6087,6170,10212,6223,6299,10583,10798)    
   descrip=' (p p-bar)->V h'//j//'-> (b b-bar) +missing Et where h'//j//' is SM-like  '       //label 
  case(10596)    
   descrip=' (p p-bar)->V h'//j//'-> (b b-bar) l nu where h'//j//' is SM-like  '       //label 
  case(6091,1268)    
   descrip=' (p p-bar)->V h'//j//'-> ll + X where h'//j//' is SM-like '       //label
  case(10010)    
   descrip=' (p p-bar)->V (h'//j//')/VBF-> (b b-bar) q q where h'//j//' is SM-like '       //label 
  case(10607)    
   descrip=' (p p-bar)->V (h'//j//')/VBF-> (b b-bar)+... where h'//j//' is SM-like '       //label 
  case(6436)    
   descrip=' (p p-bar)->V (h'//j//')-> (b b-bar)+...'       //label 
  case(9248,10133,10439,6305,6286)    
   descrip=' (p p-bar)->h'//j//'+...->tau tau +... where h'//j//' is SM-like  ' //label
  case(4800,5845,6171)    
   descrip=' (p p-bar)->h'//j//'+...->tau tau (2 jets) where h'//j//' is SM-like  ' //label  
  case(5871)  
   descrip=' (p p-bar)->h'//j//'+...->W W +... ->l l nu nu +... where h'//j//' is SM-like ' //label
  case(6082)  
   descrip=' (p p-bar)->h'//j//'+...->V V +... ->e mu missing Et +... where h'//j//' is SM-like ' //label
  case(6182,6219)  
   descrip=' (p p-bar)->h'//j//'+...->V V +... ->l l missing Et +... where h'//j//' is SM-like ' //label
  case(6276)  
   descrip=' (p p-bar)->h'//j//'+...->V V +... ->l l l missing Et +... where h'//j//' is SM-like ' //label   
  case(10573)  
   descrip=' (p p-bar)->h'//j//'+...->V V +... ->l l l l +... where h'//j//' is SM-like ' //label
  case(5984,9714,6006,9022,9023,0710,9887,4162,10102,4481,4468,6095,10432,6179,6302,10599)  
   descrip=' (p p-bar)->h'//j//'+...->W W +... where h'//j//' is SM-like ' //label 
  case(0024,5985,0968,5974,6083,4885)  
   descrip=' (p p-bar)->h'//j//'(b/b-bar)->(tau tau) (b/b-bar)   '          //label
  case(5739,10574)  
   descrip=' (p p-bar)->t t-bar h'//j//'->t t-bar b b-bar   '          //label  
  case(2012135,12025)  
   descrip=' (p p)->t t-bar h'//j//'->t t-bar b b-bar   '          //label  
  case(0611)  
   descrip=' (p p-bar)->h'//j//'->Z gamma   '          //label  
  case(10500)
   descrip=' (p p-bar)->V h'//j//'-> V tau tau ' //label
  case(1811)
   descrip=' t->(H'//j//'+)b->(2 quarks) b   '          //label 
  case(1812,2011138,2011151,2760,2013090,8353,7712,11002,11008)
   descrip=' t->(H'//j//'+)b->tau nu b   '          //label 
  case(1269,1270,2011094)
   descrip=' t->(H'//j//'+)b->(c s) b'          //label 
  case(11006,11017,110271,110272,14161,14162,5064,2012017,2011150) 
   descrip=' (p p)->h'//j//'/VBF->Z Z-> l l q q where h'//j//' is SM-like ' //label
  case(2011048,11004,11015,2011131) 
   descrip=' (p p)->h'//j//'/VBF->Z Z-> l l l l where h'//j//' is SM-like ' //label
  case(2011162,1415,2012092) 
   descrip=' (p p)->h'//j//'/VBF/V h'//j//'->Z Z-> l l l l where h'//j//' is SM-like ' //label
  case(11025,1997,12041) 
   descrip=' (p p)->h'//j//'/VBF/V/tt h'//j//'->Z Z-> l l l l where h'//j//' is SM-like ' //label
  case(20130131) 
   descrip=' (p p)->h'//j//'->Z Z-> l l l l where h'//j//' is SM-like ' //label
  case(20130132) 
   descrip=' (p p)->h'//j//'/ggF h->Z Z-> l l l l ' //label
  case(20130133) 
   descrip=' (p p)->h'//j//'/VBF/V h->Z Z-> l l l l ' //label
  case(130021) 
   descrip=' (p p)->h'//j//'->Z Z-> l l l l (low mass) where h'//j//' is SM-like ' //label
  case(130022) 
   descrip=' (p p)->h'//j//'->Z Z-> l l l l (high mass) where h'//j//' is SM-like ' //label
  case(11005,11016,11026,3478)
   descrip=' (p p)->h'//j//'/VBF->V V-> l l nu nu where h'//j//' is SM-like ' //label
  case(3357,2011148,2012016) 
   descrip=' (p p)->h'//j//'->V V-> l l nu nu where h'//j//' is SM-like ' //label
  case(11013,11028)
   descrip=' (p p)->h'//j//'/VBF->V V-> l l tau tau where h'//j//' is SM-like ' //label
  case(2011026) 
   descrip=' (p p)->h'//j//'/VBF->V V where h'//j//' is SM-like ' //label
  case(5429,2011052,2011111,2011134) 
   descrip=' (p p)->h'//j//'->W W ' //label
  case(11034,12039,13009,2012078)  
   descrip=' (p p)->W(h'//j//')->W W W where h'//j//' is SM-like ' //label    
  case(12006)  
   descrip=' (p p)->W(h'//j//')->W tau tau   ' //label
  case(12051)  
   descrip=' (p p)->V(h'//j//')->V tau tau   ' //label
  case(2012012,2012158,2013030) 
   descrip=' (p p)->h'//j//'->W W where h'//j//' is SM-like ' //label
  case(110212)
   descrip=' (p p)->V h'//j//'/VBF->gamma gamma    '    //label
  case(2011025,2011085,2011161,5895,1414,2012091,2012168,1487,12001,12015,13001,11010,11030,11021) 
   descrip=' (p p)->h'//j//'+...->gamma gamma+... where h'//j//' is SM-like ' //label
  case(13006, 13075515,2013009) 
   descrip=' (p p)->h'//j//'+...->gamma Z+... where h'//j//' is SM-like ' //label
  case(11031,12044,13012)
   descrip=' (p p)->V h'//j//'->b b where h'//j//' is SM-like ' //label
  case(13011)
   descrip=' (p p)->h'//j//'/VBF->bb+... where h'//j//' is SM-like ' //label
  case(2011020,2011021) 
   descrip=' (p p)->h'//j//'->mu mu (lower mass range)  ' //label
  case(2011005,3615,2012018) 
   descrip=' (p p)->h'//j//'/VBF->W W where h'//j//' is SM-like ' //label
  case(12046) 
   descrip=' (p p)->h'//j//'->W W-> l nu q q where h'//j//' is SM-like ' //label
  case(11003,11014,2577,11024,1489,12042,13003)  
   descrip=' (p p)->h'//j//'+...->W W +... where h'//j//' is SM-like ' //label 
  case(2748,1408,11011,2011163,11022,11032,1488,12008,12045,2011157,2011112,2011135,2012019) 
   descrip=' (p p)->h'//j//'+... where h'//j//' is SM-like ' //label        
  case(7214) 
   descrip=' (p p)->h'//j//'+... where h'//j//' is SM-like ' //label        
  case(10002,5003,2011132,2012094,110201,110292,12050) 
   descrip=' (p p)->h'//j//'->tau tau ' //label
  case(11009,11020,2011133,2012014) 
   descrip=' (p p)->h'//j//'/VBF->tau tau +... where h'//j//' is SM-like ' //label
  case(2012160,12043)
   descrip=' (p p)->h'//j//'->tau tau +... where h'//j//' is SM-like ' //label
  case(2013010)
   descrip=' (p p)->h'//j//'->mu mu +... where h'//j//' is SM-like ' //label
  case(2012015)    
   descrip=' (p p)->V h'//j//'-> (b b-bar) + X where h'//j//' is SM-like  '       //label 
  case(110291) 
   descrip=' (p p)->h'//j//'/VBF/V h'//j//'/tt h'//j//'->tau tau +... where h'//j//' is SM-like ' //label
  case(2011103,2012161,11012)  
   descrip=' (p p)->V(h'//j//')->V (b b-bar)   '  //label        
  case(13022)  
   descrip=' (p p)->h'//j//'(VBF)->WW   '  //label        
  case(13013)  
   descrip=' (p p)->h'//j//'(VBF)->V (invisible)   '  //label        
  case(2013011,13018)  
   descrip=' (p p)->Vh'//j//'->V (invisible)   '  //label        
  case default
   stop 'wrong input to function outputproc_t1 in module S95tables (1)'  
  end select   

! New description string based on data file input
! Added by OS 2012-03-12
     if(S95_t1(tlistn)%desc.NE.'') then
        descrip = trim(S95_t1(tlistn)%desc) // ', h='//j
        if (S95_t1(tlistn)%SMlike.EQ.1) then
         descrip = trim(descrip)//' where h is SM-like'
        endif
         descrip = trim(descrip)//' '//label
     endif
 end subroutine outputproc_t1
 
 !********************************************************      
 subroutine outputproc_t2(tlistn,ii,jj,k,descrip)
 !********************************************************            
 ! uses information about the process to output a description
 ! for processes using table type 1
 ! note: at the moment, np(x) (and so ii and jj) needs to be 1 digit long i.e. np(x)<10
 !******************************************************** 
  implicit none
  !--------------------------------------input
  integer :: tlistn
  integer :: ii,jj,k
  !-----------------------------------internal
  character(LEN=1) :: j,i
  character(LEN=45) :: label
  character(LEN=200):: descrip
  !-------------------------------------------
               
  if((ii.ne.0).and.(jj.ne.0))then      
   write(i,'(I1)')ii 
   write(j,'(I1)')jj      
  else
   i='i'
   j='j'
  endif   

  if(k.eq.21)then
   label=''  !no need to lable each line in Key.dat
  else
   label='('//trim(S95_t2(tlistn)%label)//')'
  endif
                         
  select case(S95_t2(tlistn)%id)
  case(150)      
   descrip=' (ee)->(h'//j//'->h'//i//' h'//i//')Z->(b b b b)Z   '         //label      
  case(160)      
   descrip=' (ee)->(h'//j//'->h'//i//' h'//i//')Z->(tau tau tau tau)Z   ' //label     
  case(180)      
   descrip=' (ee)->(h'//j//' h'//i//')->(b b b b)   '                     //label
  case(190)      
   descrip=' (ee)->(h'//j//' h'//i//')->(tau tau tau tau)   '             //label
  case(200)      
   descrip=' (ee)->(h'//j//'->h'//i//' h'//i//')h'//i//'->(b b b b)b b   '//label
  case(210)      
   descrip=' (ee)->(h'//j//'->h'//i//' h'//i//')h'//i//'->(tau tau tau tau)tau tau   '//label
  case(220)      
   descrip=' (ee)->(h'//j//'->h'//i//' h'//i//')Z->(b b)(tau tau)Z   '    //label
  case(230)      
   descrip=' (ee)->(h'//j//'->b b)(h'//i//'->tau tau)   '                 //label
  case(240)      
   descrip=' (ee)->(h'//j//'->tau tau)(h'//i//'->b b)   '                 //label
  case(905)
   descrip=' (ee)->(C'//j//'+)(C'//j//'-)-> (q q N'//i//') (q q N'//i//')   ' //label
  case(906)
   descrip=' (ee)->(C'//j//'+)(C'//j//'-)-> q q l nu N'//i//' N'//i//'   ' //label
  case(907)
   descrip=' (ee)->(C'//j//'+)(C'//j//'-)-> (l nu N'//i//') (l nu N'//i//')   ' //label
  case(908)
   descrip=' (ee)->(C'//j//'+)(C'//j//'-) with all C'//j//' decaying to W + N'//i//' ' //label
  case(909)
   descrip=' (ee)->(N'//j//') N'//i//'-> (q q N'//i//') N'//i//'   ' //label
  case(910)      
   descrip=' (ee)->N'//j//' N'//i//' with all N'//j//' decaying to Z + N'//i//'  ' //label
  case(3381)  
   descrip=' (p p-bar)->h'//j//'->h'//i//' h'//i//'->mu mu mu mu   '          //label  
  case(3382)  
   descrip=' (p p-bar)->h'//j//'->h'//i//' h'//i//'->tau tau mu mu   '          //label
  case(6227)  
   descrip=' (p p-bar)->h'//j//'(b/b-bar)->(b b-bar) (b/b-bar) or (tau tau) (b/b-bar) '     //label 
  case default
   stop 'wrong input to function outputproc_t2 in module S95tables (2)' 
  end select        
            
 end subroutine outputproc_t2
 
 !******************************************************************
 subroutine model_likeness(j,id,t,model_like,sigmaXbr)
 !***************************************************************** 
 ! Tests how Standard Model-like a parameter point is            
 ! 0 means Mi.ge.MSingleLim (treat as single channel)
 ! 1 passes the SM-like test and Mi.lt.MSingleLim
 ! -1 fails the SM-like test and Mi.lt.MSingleLim 
  use usefulbits, only : dataset,div, vsmall, iselementofarray
  use theory_BRfunctions
  use theory_XS_SM_functions   
  implicit none      
  !--------------------------------------input      
  type(dataset) :: t 
  integer :: id,j
  !-----------------------------------internal 
  integer :: ns,nb,n
  !--TS 14/03/2011: For revamped model-likeness test method
  double precision,allocatable :: channel_rat(:,:), channel_SM(:,:)
  double precision,allocatable :: XS_SM_temp(:), BR_SM_temp(:)
  double precision :: SMrate, weight, c
  integer :: ic,nc,nc_rel
  !----
  double precision,allocatable :: XS_rat(:), BR_rat(:)
  integer :: model_like,testSMratios
  double precision :: sigmaXbr 
  integer :: is,ib
  double precision :: s,b
  double precision,allocatable :: dsbys(:),dbbyb(:),dcbyc(:)
  logical :: correct_properties
  double precision,parameter :: unset=-9.9999D6

  correct_properties=.True.
  ns=-1
  nb=-1
  nc=-1

  n=t1elementnumberfromid(S95_t1,id) 
  select case(id)
  case(711,713,721,723,731,733,741,743,5739,10574,6224,6225,6226,6276,6301,6309)
   !these have a very simple model-likeness test, so we can have a non-zero deltax
  case default
   if(S95_t1(n)%deltax.gt.0.0D0)then
     write(*,*)'hello id=',id,'deltax=',S95_t1(n)%deltax
     stop'error in subroutine model_likeness (1)'
   endif
  end select

  select case(id)
  case(8961,0598)
!   ns = 3; nb = 2; call initialise_XS_rat_BR_rat   
   nc = 6; call initialise_channel_rat_SM
        
!   XS_rat(1) = t%tev%XS_hjW_ratio(j)
!   XS_rat(2) = t%tev%XS_hj_ratio(j)  
!   XS_rat(3) = t%tev%XS_hjZ_ratio(j) 
  
!   BR_rat(1) = div(t%BR_hjbb(j) , t%BR_Hbb_SM(j) ,0.0D0,1.0D0) 
!   BR_rat(2) = div(t%BR_hjWW(j) , t%BR_HWW_SM(j) ,0.0D0,1.0D0)     

   channel_rat(1,:) = (/ t%tev%XS_hj_ratio(j)  , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(2,:) = (/ t%tev%XS_hjW_ratio(j) , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(3,:) = (/ t%tev%XS_hjZ_ratio(j) , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(4,:) = (/ t%tev%XS_hj_ratio(j)  , div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0) /)
   channel_rat(5,:) = (/ t%tev%XS_hjW_ratio(j) , div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0) /)
   channel_rat(6,:) = (/ t%tev%XS_hjZ_ratio(j) , div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0) /)

   channel_SM(1,:) = (/ t%tev%XS_H_SM(j)   , t%BR_HWW_SM(j) /)
   channel_SM(2,:) = (/ t%tev%XS_HW_SM(j) , t%BR_HWW_SM(j) /) 
   channel_SM(3,:) = (/ t%tev%XS_HZ_SM(j) , t%BR_HWW_SM(j) /)
   channel_SM(4,:) = (/ t%tev%XS_H_SM(j)   , t%BR_Hbb_SM(j) /)
   channel_SM(5,:) = (/ t%tev%XS_HW_SM(j) , t%BR_Hbb_SM(j) /) 
   channel_SM(6,:) = (/ t%tev%XS_HZ_SM(j) , t%BR_Hbb_SM(j) /)  
      
  case(10010,10607)
!   ns = 4; nb = 2; call initialise_XS_rat_BR_rat      
   nc = 8; call initialise_channel_rat_SM
        
!   XS_rat(1) = t%tev%XS_hjW_ratio(j)
!   XS_rat(2) = t%tev%XS_vbf_ratio(j)  
!   XS_rat(3) = t%tev%XS_hjZ_ratio(j) 
!   XS_rat(4) = t%tev%XS_tthj_ratio(j)

!   BR_rat(1) = div(t%BR_hjbb(j) , t%BR_Hbb_SM(j),0.0D0,1.0D0)
!   BR_rat(2) = div(( t%BR_hjss(j)+t%BR_hjcc(j)+t%BR_hjbb(j)+t%BR_hjgg(j) ) &
!             & , t%BR_Hjets_SM(j),0.0D0,1.0D0)


   channel_rat(1,:) = (/ t%tev%XS_hj_ratio(j)  , &
   & div((t%BR_hjss(j)+t%BR_hjcc(j)+t%BR_hjbb(j)+t%BR_hjgg(j)),t%BR_Hjets_SM(j),0.0D0,1.0D0) /)
   channel_rat(2,:) = (/ t%tev%XS_vbf_ratio(j) , &
   & div((t%BR_hjss(j)+t%BR_hjcc(j)+t%BR_hjbb(j)+t%BR_hjgg(j)),t%BR_Hjets_SM(j),0.0D0,1.0D0) /)
   channel_rat(3,:) = (/ t%tev%XS_hjW_ratio(j) , &
   & div((t%BR_hjss(j)+t%BR_hjcc(j)+t%BR_hjbb(j)+t%BR_hjgg(j)),t%BR_Hjets_SM(j),0.0D0,1.0D0) /)
   channel_rat(4,:) = (/ t%tev%XS_hjZ_ratio(j) , &
   & div((t%BR_hjss(j)+t%BR_hjcc(j)+t%BR_hjbb(j)+t%BR_hjgg(j)),t%BR_Hjets_SM(j),0.0D0,1.0D0) /)
   channel_rat(5,:) = (/ t%tev%XS_hj_ratio(j)  , div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0) /)
   channel_rat(6,:) = (/ t%tev%XS_vbf_ratio(j) , div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0) /)
   channel_rat(7,:) = (/ t%tev%XS_hjW_ratio(j) , div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0) /)
   channel_rat(8,:) = (/ t%tev%XS_hjZ_ratio(j) , div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0) /)

   channel_SM(1,:) = (/ t%tev%XS_H_SM(j)   , t%BR_Hjets_SM(j) /)
   channel_SM(2,:) = (/ t%tev%XS_vbf_SM(j) , t%BR_Hjets_SM(j) /)
   channel_SM(3,:) = (/ t%tev%XS_HW_SM(j) , t%BR_Hjets_SM(j) /) 
   channel_SM(4,:) = (/ t%tev%XS_HZ_SM(j) , t%BR_Hjets_SM(j) /)
   channel_SM(5,:) = (/ t%tev%XS_H_SM(j)   , t%BR_Hbb_SM(j) /)
   channel_SM(6,:) = (/ t%tev%XS_vbf_SM(j) , t%BR_Hbb_SM(j) /)
   channel_SM(7,:) = (/ t%tev%XS_HW_SM(j) , t%BR_Hbb_SM(j) /) 
   channel_SM(8,:) = (/ t%tev%XS_HZ_SM(j) , t%BR_Hbb_SM(j) /)  

 case(6436)
   nc = 2; call initialise_channel_rat_SM
        
   channel_rat(1,:) = (/ t%tev%XS_hjW_ratio(j) , div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0) /)
   channel_rat(2,:) = (/ t%tev%XS_hjZ_ratio(j) , div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0) /)

   channel_SM(1,:) = (/ t%tev%XS_HW_SM(j) , t%BR_Hbb_SM(j) /) 
   channel_SM(2,:) = (/ t%tev%XS_HZ_SM(j) , t%BR_Hbb_SM(j) /)  

  case(9290)
!   ns = 4; nb = 4; call initialise_XS_rat_BR_rat     
   nc = 16; call initialise_channel_rat_SM
        
!   XS_rat(1) = t%tev%XS_hjW_ratio(j)
!   XS_rat(2) = t%tev%XS_hj_ratio(j)  
!   XS_rat(3) = t%tev%XS_hjZ_ratio(j) 
!   XS_rat(4) = t%tev%XS_vbf_ratio(j) 
  
!   BR_rat(1) = div(    t%BR_hjbb(j) , t%BR_Hbb_SM(j)    ,0.0D0,1.0D0)
!   BR_rat(2) = div(    t%BR_hjWW(j) , t%BR_HWW_SM(j)    ,0.0D0,1.0D0)   
!   BR_rat(3) = div(t%BR_hjtautau(j) , t%BR_Htautau_SM(j),0.0D0,1.0D0) 
!   BR_rat(4) = div(  t%BR_hjgaga(j) , t%BR_Hgaga_SM(j)  ,0.0D0,1.0D0)

   channel_rat(1,:) = (/ t%tev%XS_hj_ratio(j)  , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(2,:) = (/ t%tev%XS_vbf_ratio(j) , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(3,:) = (/ t%tev%XS_hjW_ratio(j) , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(4,:) = (/ t%tev%XS_hjZ_ratio(j) , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(5,:) = (/ t%tev%XS_hj_ratio(j)  , div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(6,:) = (/ t%tev%XS_vbf_ratio(j) , div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(7,:) = (/ t%tev%XS_hjW_ratio(j) , div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(8,:) = (/ t%tev%XS_hjZ_ratio(j) , div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(9,:) = (/ t%tev%XS_hj_ratio(j)  , div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0) /)
   channel_rat(10,:) = (/ t%tev%XS_vbf_ratio(j) , div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0) /)
   channel_rat(11,:) = (/ t%tev%XS_hjW_ratio(j) , div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0) /)
   channel_rat(12,:) = (/ t%tev%XS_hjZ_ratio(j) , div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0) /)
   channel_rat(13,:) = (/ t%tev%XS_hj_ratio(j)  , div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)
   channel_rat(14,:) = (/ t%tev%XS_vbf_ratio(j) , div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)
   channel_rat(15,:) = (/ t%tev%XS_hjW_ratio(j) , div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)
   channel_rat(16,:) = (/ t%tev%XS_hjZ_ratio(j) , div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)
   
   channel_SM(1,:) = (/ t%tev%XS_H_SM(j)   , t%BR_HWW_SM(j) /)
   channel_SM(2,:) = (/ t%tev%XS_vbf_SM(j) , t%BR_HWW_SM(j) /)
   channel_SM(3,:) = (/ t%tev%XS_HW_SM(j) , t%BR_HWW_SM(j) /) 
   channel_SM(4,:) = (/ t%tev%XS_HZ_SM(j) , t%BR_HWW_SM(j) /)
   channel_SM(5,:) = (/ t%tev%XS_H_SM(j)   , t%BR_Htautau_SM(j) /)
   channel_SM(6,:) = (/ t%tev%XS_vbf_SM(j) , t%BR_Htautau_SM(j) /)
   channel_SM(7,:) = (/ t%tev%XS_HW_SM(j) , t%BR_Htautau_SM(j) /) 
   channel_SM(8,:) = (/ t%tev%XS_HZ_SM(j) , t%BR_Htautau_SM(j) /)
   channel_SM(9,:) = (/ t%tev%XS_H_SM(j)   , t%BR_Hbb_SM(j) /)
   channel_SM(10,:) = (/ t%tev%XS_vbf_SM(j) , t%BR_Hbb_SM(j) /)
   channel_SM(11,:) = (/ t%tev%XS_HW_SM(j) , t%BR_Hbb_SM(j) /) 
   channel_SM(12,:) = (/ t%tev%XS_HZ_SM(j) , t%BR_Hbb_SM(j) /)  
   channel_SM(13,:) = (/ t%tev%XS_H_SM(j)   , t%BR_Hgaga_SM(j) /)
   channel_SM(14,:) = (/ t%tev%XS_vbf_SM(j) , t%BR_Hgaga_SM(j) /)
   channel_SM(15,:) = (/ t%tev%XS_HW_SM(j) , t%BR_Hgaga_SM(j) /) 
   channel_SM(16,:) = (/ t%tev%XS_HZ_SM(j) , t%BR_Hgaga_SM(j) /)  

     
  case(9674,9897,9999)
!   ns = 4; nb = 3; call initialise_XS_rat_BR_rat  
   nc = 12; call initialise_channel_rat_SM
        
!   XS_rat(1) = t%tev%XS_hjW_ratio(j)
!   XS_rat(2) = t%tev%XS_hj_ratio(j)  
!   XS_rat(3) = t%tev%XS_hjZ_ratio(j) 
!   XS_rat(4) = t%tev%XS_vbf_ratio(j) 
  
!   BR_rat(1) = div(    t%BR_hjbb(j) , t%BR_Hbb_SM(j)     ,0.0D0,1.0D0)
!   BR_rat(2) = div(    t%BR_hjWW(j) , t%BR_HWW_SM(j)     ,0.0D0,1.0D0)   
!   BR_rat(3) = div(t%BR_hjtautau(j) , t%BR_Htautau_SM(j) ,0.0D0,1.0D0) 

   channel_rat(1,:) = (/ t%tev%XS_hj_ratio(j)  , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(2,:) = (/ t%tev%XS_vbf_ratio(j) , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(3,:) = (/ t%tev%XS_hjW_ratio(j) , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(4,:) = (/ t%tev%XS_hjZ_ratio(j) , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(5,:) = (/ t%tev%XS_hj_ratio(j)  , div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(6,:) = (/ t%tev%XS_vbf_ratio(j) , div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(7,:) = (/ t%tev%XS_hjW_ratio(j) , div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(8,:) = (/ t%tev%XS_hjZ_ratio(j) , div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(9,:) = (/ t%tev%XS_hj_ratio(j)  , div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0) /)
   channel_rat(10,:) = (/ t%tev%XS_vbf_ratio(j) , div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0) /)
   channel_rat(11,:) = (/ t%tev%XS_hjW_ratio(j) , div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0) /)
   channel_rat(12,:) = (/ t%tev%XS_hjZ_ratio(j) , div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0) /)
   
   channel_SM(1,:) = (/ t%tev%XS_H_SM(j)   , t%BR_HWW_SM(j) /)
   channel_SM(2,:) = (/ t%tev%XS_vbf_SM(j) , t%BR_HWW_SM(j) /)
   channel_SM(3,:) = (/ t%tev%XS_HW_SM(j) , t%BR_HWW_SM(j) /) 
   channel_SM(4,:) = (/ t%tev%XS_HZ_SM(j) , t%BR_HWW_SM(j) /)
   channel_SM(5,:) = (/ t%tev%XS_H_SM(j)   , t%BR_Htautau_SM(j) /)
   channel_SM(6,:) = (/ t%tev%XS_vbf_SM(j) , t%BR_Htautau_SM(j) /)
   channel_SM(7,:) = (/ t%tev%XS_HW_SM(j) , t%BR_Htautau_SM(j) /) 
   channel_SM(8,:) = (/ t%tev%XS_HZ_SM(j) , t%BR_Htautau_SM(j) /)
   channel_SM(9,:) = (/ t%tev%XS_H_SM(j)   , t%BR_Hbb_SM(j) /)
   channel_SM(10,:) = (/ t%tev%XS_vbf_SM(j) , t%BR_Hbb_SM(j) /)
   channel_SM(11,:) = (/ t%tev%XS_HW_SM(j) , t%BR_Hbb_SM(j) /) 
   channel_SM(12,:) = (/ t%tev%XS_HZ_SM(j) , t%BR_Hbb_SM(j) /)   
     
  case(7081,9166,9483,5586,9642,1266,0432,9891,5285,3935,6087,6170,10212,6223,6299,10583,10798,10596)
!   ns = 2; nb = 1; call initialise_XS_rat_BR_rat      
   nc = 2; call initialise_channel_rat_SM
  
!   XS_rat(1) = t%tev%XS_hjW_ratio(j)
!   XS_rat(2) = t%tev%XS_hjZ_ratio(j)  
  
!   BR_rat(1) = div(t%BR_hjbb(j) , t%BR_Hbb_SM(j),0.0D0,1.0D0)
 
   channel_rat(1,:) = (/ t%tev%XS_hjW_ratio(j) , div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0) /)
   channel_rat(2,:) = (/ t%tev%XS_hjZ_ratio(j) , div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0) /) 

   channel_SM(1,:) = (/ t%tev%XS_HW_SM(j) , t%BR_Hbb_SM(j) /) 
   channel_SM(2,:) = (/ t%tev%XS_HZ_SM(j) , t%BR_Hbb_SM(j) /)  
   

     
  case(10500)
!   ns = 2; nb = 1; call initialise_XS_rat_BR_rat      
   nc = 2; call initialise_channel_rat_SM
        
!   XS_rat(1) = t%tev%XS_hjW_ratio(j)
!   XS_rat(2) = t%tev%XS_hjZ_ratio(j)  
  
!   BR_rat(1) = div(t%BR_hjtautau(j) , t%BR_Htautau_SM(j) ,0.0D0,1.0D0) 

   channel_rat(1,:) = (/ t%tev%XS_hjW_ratio(j) , div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(2,:) = (/ t%tev%XS_hjZ_ratio(j) , div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /) 

   channel_SM(1,:) = (/ t%tev%XS_HW_SM(j) , t%BR_Htautau_SM(j) /) 
   channel_SM(2,:) = (/ t%tev%XS_HZ_SM(j) , t%BR_Htautau_SM(j) /)

  case(9248,10133,10439)
!   ns = 4; nb = 1; call initialise_XS_rat_BR_rat    
   nc = 4; call initialise_channel_rat_SM

!   XS_rat(1) = t%tev%XS_hjW_ratio(j)
!   XS_rat(2) = t%tev%XS_hj_ratio(j)  
!   XS_rat(3) = t%tev%XS_hjZ_ratio(j)  
!   XS_rat(4) = t%tev%XS_vbf_ratio(j)
  
!   BR_rat(1) = div(t%BR_hjtautau(j) , t%BR_Htautau_SM(j),0.0D0,1.0D0)

   channel_rat(1,:) = (/ t%tev%XS_hj_ratio(j)  , div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(2,:) = (/ t%tev%XS_vbf_ratio(j) , div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(3,:) = (/ t%tev%XS_hjW_ratio(j) , div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(4,:) = (/ t%tev%XS_hjZ_ratio(j) , div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /) 
   
   channel_SM(1,:) = (/ t%tev%XS_H_SM(j)   , t%BR_Htautau_SM(j) /)
   channel_SM(2,:) = (/ t%tev%XS_vbf_SM(j) , t%BR_Htautau_SM(j) /)
   channel_SM(3,:) = (/ t%tev%XS_HW_SM(j) , t%BR_Htautau_SM(j) /) 
   channel_SM(4,:) = (/ t%tev%XS_HZ_SM(j) , t%BR_Htautau_SM(j) /)

  case(4800)
!   ns = 4; nb = 3; call initialise_XS_rat_BR_rat
   nc = 12; call initialise_channel_rat_SM
        
!   XS_rat(1) = t%tev%XS_hjW_ratio(j)
!   XS_rat(2) = t%tev%XS_hj_ratio(j)  
!   XS_rat(3) = t%tev%XS_hjZ_ratio(j)  
!   XS_rat(4) = t%tev%XS_vbf_ratio(j)
  
!   BR_rat(1) = div(t%BR_hjtautau(j) , t%BR_Htautau_SM(j) ,0.0D0,1.0D0) 
!   BR_rat(2) = div(t%BR_hjbb(j) , t%BR_Hbb_SM(j),0.0D0,1.0D0)
!   BR_rat(3) = div(( t%BR_hjss(j)+t%BR_hjcc(j)+t%BR_hjbb(j)+t%BR_hjgg(j) ) &
!             & , t%BR_Hjets_SM(j),0.0D0,1.0D0)

   channel_rat(1,:) = (/ t%tev%XS_hj_ratio(j)  , div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(2,:) = (/ t%tev%XS_vbf_ratio(j) , div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(3,:) = (/ t%tev%XS_hjW_ratio(j) , div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(4,:) = (/ t%tev%XS_hjZ_ratio(j) , div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(5,:) = (/ t%tev%XS_hj_ratio(j)  , &
   & div((t%BR_hjss(j)+t%BR_hjcc(j)+t%BR_hjbb(j)+t%BR_hjgg(j)),t%BR_Hjets_SM(j),0.0D0,1.0D0) /)
   channel_rat(6,:) = (/ t%tev%XS_vbf_ratio(j) , &
   & div((t%BR_hjss(j)+t%BR_hjcc(j)+t%BR_hjbb(j)+t%BR_hjgg(j)),t%BR_Hjets_SM(j),0.0D0,1.0D0) /)
   channel_rat(7,:) = (/ t%tev%XS_hjW_ratio(j) , &
   & div((t%BR_hjss(j)+t%BR_hjcc(j)+t%BR_hjbb(j)+t%BR_hjgg(j)),t%BR_Hjets_SM(j),0.0D0,1.0D0) /)
   channel_rat(8,:) = (/ t%tev%XS_hjZ_ratio(j) , &
   & div((t%BR_hjss(j)+t%BR_hjcc(j)+t%BR_hjbb(j)+t%BR_hjgg(j)),t%BR_Hjets_SM(j),0.0D0,1.0D0) /)
   channel_rat(9,:) = (/ t%tev%XS_hj_ratio(j)  , div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0) /)
   channel_rat(10,:) = (/ t%tev%XS_vbf_ratio(j) , div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0) /)
   channel_rat(11,:) = (/ t%tev%XS_hjW_ratio(j) , div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0) /)
   channel_rat(12,:) = (/ t%tev%XS_hjZ_ratio(j) , div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0) /)

   channel_SM(1,:) = (/ t%tev%XS_H_SM(j)   , t%BR_Htautau_SM(j) /)
   channel_SM(2,:) = (/ t%tev%XS_vbf_SM(j) , t%BR_Htautau_SM(j) /)
   channel_SM(3,:) = (/ t%tev%XS_HW_SM(j) , t%BR_Htautau_SM(j) /) 
   channel_SM(4,:) = (/ t%tev%XS_HZ_SM(j) , t%BR_Htautau_SM(j) /)
   channel_SM(5,:) = (/ t%tev%XS_H_SM(j)   , t%BR_Hjets_SM(j) /)
   channel_SM(6,:) = (/ t%tev%XS_vbf_SM(j) , t%BR_Hjets_SM(j) /)
   channel_SM(7,:) = (/ t%tev%XS_HW_SM(j) , t%BR_Hjets_SM(j) /) 
   channel_SM(8,:) = (/ t%tev%XS_HZ_SM(j) , t%BR_Hjets_SM(j) /)
   channel_SM(9,:) = (/ t%tev%XS_H_SM(j)   , t%BR_Hbb_SM(j) /)
   channel_SM(10,:) = (/ t%tev%XS_vbf_SM(j) , t%BR_Hbb_SM(j) /)
   channel_SM(11,:) = (/ t%tev%XS_HW_SM(j) , t%BR_Hbb_SM(j) /) 
   channel_SM(12,:) = (/ t%tev%XS_HZ_SM(j) , t%BR_Hbb_SM(j) /)       

  case(5845)
!   ns = 4; nb = 2; call initialise_XS_rat_BR_rat 
   nc = 8; call initialise_channel_rat_SM
        
!   XS_rat(1) = t%tev%XS_hjW_ratio(j)
!   XS_rat(2) = t%tev%XS_hj_ratio(j)  
!   XS_rat(3) = t%tev%XS_hjZ_ratio(j)  
!   XS_rat(4) = t%tev%XS_vbf_ratio(j)
  
!   BR_rat(1) = div(t%BR_hjtautau(j) , t%BR_Htautau_SM(j)   ,0.0D0,1.0D0)  
!   BR_rat(2) = div(( t%BR_hjss(j)+t%BR_hjcc(j)+t%BR_hjbb(j)+t%BR_hjgg(j) ) &
!             & , t%BR_Hjets_SM(j) ,0.0D0,1.0D0)  

   channel_rat(1,:) = (/ t%tev%XS_hj_ratio(j)  , div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(2,:) = (/ t%tev%XS_vbf_ratio(j) , div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(3,:) = (/ t%tev%XS_hjW_ratio(j) , div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(4,:) = (/ t%tev%XS_hjZ_ratio(j) , div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(5,:) = (/ t%tev%XS_hj_ratio(j)  , &
   & div((t%BR_hjss(j)+t%BR_hjcc(j)+t%BR_hjbb(j)+t%BR_hjgg(j)),t%BR_Hjets_SM(j),0.0D0,1.0D0) /)
   channel_rat(6,:) = (/ t%tev%XS_vbf_ratio(j) , &
   & div((t%BR_hjss(j)+t%BR_hjcc(j)+t%BR_hjbb(j)+t%BR_hjgg(j)),t%BR_Hjets_SM(j),0.0D0,1.0D0) /)
   channel_rat(7,:) = (/ t%tev%XS_hjW_ratio(j) , &
   & div((t%BR_hjss(j)+t%BR_hjcc(j)+t%BR_hjbb(j)+t%BR_hjgg(j)),t%BR_Hjets_SM(j),0.0D0,1.0D0) /)
   channel_rat(8,:) = (/ t%tev%XS_hjZ_ratio(j) , &
   & div((t%BR_hjss(j)+t%BR_hjcc(j)+t%BR_hjbb(j)+t%BR_hjgg(j)),t%BR_Hjets_SM(j),0.0D0,1.0D0) /)

   channel_SM(1,:) = (/ t%tev%XS_H_SM(j)   , t%BR_Htautau_SM(j) /)
   channel_SM(2,:) = (/ t%tev%XS_vbf_SM(j) , t%BR_Htautau_SM(j) /)
   channel_SM(3,:) = (/ t%tev%XS_HW_SM(j) , t%BR_Htautau_SM(j) /) 
   channel_SM(4,:) = (/ t%tev%XS_HZ_SM(j) , t%BR_Htautau_SM(j) /)
   channel_SM(5,:) = (/ t%tev%XS_H_SM(j)   , t%BR_Hjets_SM(j) /)
   channel_SM(6,:) = (/ t%tev%XS_vbf_SM(j) , t%BR_Hjets_SM(j) /)
   channel_SM(7,:) = (/ t%tev%XS_HW_SM(j) , t%BR_Hjets_SM(j) /) 
   channel_SM(8,:) = (/ t%tev%XS_HZ_SM(j) , t%BR_Hjets_SM(j) /)


  case(5858,6177,6295,1887,10065,10485,4960)
!   ns = 4; nb = 1; call initialise_XS_rat_BR_rat   
   nc = 4; call initialise_channel_rat_SM
        
!   XS_rat(1) = t%tev%XS_hjW_ratio(j)
!   XS_rat(2) = t%tev%XS_hj_ratio(j)  
!   XS_rat(3) = t%tev%XS_hjZ_ratio(j)  
!   XS_rat(4) = t%tev%XS_vbf_ratio(j)

!   BR_rat(1) = div(t%BR_hjgaga(j) , t%BR_Hgaga_SM(j),0.0D0,1.0D0)  

   channel_rat(1,:) = (/ t%tev%XS_hj_ratio(j)  , div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)
   channel_rat(2,:) = (/ t%tev%XS_vbf_ratio(j) , div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)
   channel_rat(3,:) = (/ t%tev%XS_hjW_ratio(j) , div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)
   channel_rat(4,:) = (/ t%tev%XS_hjZ_ratio(j) , div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /) 
   
   channel_SM(1,:) = (/ t%tev%XS_H_SM(j)   , t%BR_Hgaga_SM(j) /)
   channel_SM(2,:) = (/ t%tev%XS_vbf_SM(j) , t%BR_Hgaga_SM(j) /)
   channel_SM(3,:) = (/ t%tev%XS_HW_SM(j) , t%BR_Hgaga_SM(j) /) 
   channel_SM(4,:) = (/ t%tev%XS_HZ_SM(j) , t%BR_Hgaga_SM(j) /)


  case(9465,5871,9022,9023,0710,9887,5984,9714,4162,10102,6006,4481,4468,6095,10432,6179,6302,10599)
!   ns = 4; nb = 1; call initialise_XS_rat_BR_rat  
   nc = 4; call initialise_channel_rat_SM
           
!   XS_rat(1) = t%tev%XS_hjW_ratio(j)
!   XS_rat(2) = t%tev%XS_hj_ratio(j)  
!   XS_rat(3) = t%tev%XS_hjZ_ratio(j)  
!   XS_rat(4) = t%tev%XS_vbf_ratio(j)
  
!   BR_rat(1) = div( t%BR_hjWW(j) , t%BR_HWW_SM(j)  ,0.0D0,1.0D0)     

   channel_rat(1,:) = (/ t%tev%XS_hj_ratio(j)  , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(2,:) = (/ t%tev%XS_vbf_ratio(j) , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(3,:) = (/ t%tev%XS_hjW_ratio(j) , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(4,:) = (/ t%tev%XS_hjZ_ratio(j) , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /) 
   
   channel_SM(1,:) = (/ t%tev%XS_H_SM(j)   , t%BR_HWW_SM(j) /)
   channel_SM(2,:) = (/ t%tev%XS_vbf_SM(j) , t%BR_HWW_SM(j) /)
   channel_SM(3,:) = (/ t%tev%XS_HW_SM(j) , t%BR_HWW_SM(j) /) 
   channel_SM(4,:) = (/ t%tev%XS_HZ_SM(j) , t%BR_HWW_SM(j) /)
   
  case(6082,6182,6219,10573,6276)
!   ns = 4; nb = 2; call initialise_XS_rat_BR_rat
   nc = 8; call initialise_channel_rat_SM
        
!   XS_rat(1) = t%tev%XS_hjW_ratio(j)
!   XS_rat(2) = t%tev%XS_hj_ratio(j)  
!   XS_rat(3) = t%tev%XS_hjZ_ratio(j)  
!   XS_rat(4) = t%tev%XS_vbf_ratio(j)
  
!   BR_rat(1) = div(t%BR_hjWW(j) , t%BR_HWW_SM(j)  ,0.0D0,1.0D0)     
!   BR_rat(2) = div(t%BR_hjZZ(j) , t%BR_HZZ_SM(j)  ,0.0D0,1.0D0) 
   
   channel_rat(1,:) = (/ t%tev%XS_hj_ratio(j)  , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(2,:) = (/ t%tev%XS_vbf_ratio(j) , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(3,:) = (/ t%tev%XS_hjW_ratio(j) , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(4,:) = (/ t%tev%XS_hjZ_ratio(j) , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)   
   channel_rat(5,:) = (/ t%tev%XS_hj_ratio(j)  , div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(6,:) = (/ t%tev%XS_vbf_ratio(j) , div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(7,:) = (/ t%tev%XS_hjW_ratio(j) , div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(8,:) = (/ t%tev%XS_hjZ_ratio(j) , div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   
   channel_SM(1,:) = (/ t%tev%XS_H_SM(j)   , t%BR_HWW_SM(j) /)
   channel_SM(2,:) = (/ t%tev%XS_vbf_SM(j) , t%BR_HWW_SM(j) /)
   channel_SM(3,:) = (/ t%tev%XS_HW_SM(j) , t%BR_HWW_SM(j) /) 
   channel_SM(4,:) = (/ t%tev%XS_HZ_SM(j) , t%BR_HWW_SM(j) /)
   channel_SM(5,:) = (/ t%tev%XS_H_SM(j)   , t%BR_HZZ_SM(j) /)
   channel_SM(6,:) = (/ t%tev%XS_vbf_SM(j) , t%BR_HZZ_SM(j) /)
   channel_SM(7,:) = (/ t%tev%XS_HW_SM(j) , t%BR_HZZ_SM(j) /) 
   channel_SM(8,:) = (/ t%tev%XS_HZ_SM(j) , t%BR_HZZ_SM(j) /)   
   
  case(3331)
!   ns = 1; nb = 2; call initialise_XS_rat_BR_rat    
   nc = 2; call initialise_channel_rat_SM        
   
!   XS_rat(1) = t%tev%XS_hj_ratio(j)
      
!   BR_rat(1) = div( t%BR_hjWW(j) , t%BR_HWW_SM(j)  ,0.0D0,1.0D0) 
!   BR_rat(2) = div( t%BR_hjZZ(j) , t%BR_HZZ_SM(j)  ,0.0D0,1.0D0) 

   channel_rat(1,:) = (/ t%tev%XS_hj_ratio(j)  , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(2,:) = (/ t%tev%XS_hj_ratio(j)  , div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)

   channel_SM(1,:) = (/ t%tev%XS_H_SM(j)   , t%BR_HWW_SM(j) /)
   channel_SM(2,:) = (/ t%tev%XS_H_SM(j)   , t%BR_HZZ_SM(j) /)

  case(6301)
!   ns = 2; nb = 1; call initialise_XS_rat_BR_rat    
   nc = 2; call initialise_channel_rat_SM
        
!   XS_rat(1) = t%tev%XS_hjW_ratio(j)
!   XS_rat(2) = t%tev%XS_hjZ_ratio(j) 
      
!   BR_rat(1) = div( t%BR_hjWW(j) , t%BR_HWW_SM(j)  ,0.0D0,1.0D0) 

   channel_rat(1,:) = (/ t%tev%XS_hjW_ratio(j) , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(2,:) = (/ t%tev%XS_hjZ_ratio(j) , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)

   channel_SM(1,:) = (/ t%tev%XS_HW_SM(j) , t%BR_HWW_SM(j) /) 
   channel_SM(2,:) = (/ t%tev%XS_HZ_SM(j) , t%BR_HWW_SM(j) /)

  case(6309)
!   ns = 4; nb = 2; call initialise_XS_rat_BR_rat    
   nc = 5; call initialise_channel_rat_SM
        
   channel_rat(1,:) = (/ t%tev%XS_hjW_ratio(j) , div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0) /)
   channel_rat(2,:) = (/ t%tev%XS_hj_ratio(j)  , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(3,:) = (/ t%tev%XS_vbf_ratio(j) , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)   
   channel_rat(4,:) = (/ t%tev%XS_hjW_ratio(j) , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)   
   channel_rat(5,:) = (/ t%tev%XS_hjZ_ratio(j) , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)

   channel_SM(1,:) = (/ t%tev%XS_HW_SM(j) , t%BR_Hbb_SM(j) /) 
   channel_SM(2,:) = (/ t%tev%XS_H_SM(j)   , t%BR_HWW_SM(j) /)
   channel_SM(3,:) = (/ t%tev%XS_vbf_SM(j) , t%BR_HWW_SM(j) /)
   channel_SM(4,:) = (/ t%tev%XS_HW_SM(j) , t%BR_HWW_SM(j) /) 
   channel_SM(5,:) = (/ t%tev%XS_HZ_SM(j) , t%BR_HWW_SM(j) /)
  
  case(1268,6091)
!   ns = 2; nb = 2; call initialise_XS_rat_BR_rat    
   nc = 4; call initialise_channel_rat_SM
        
!   XS_rat(1) = t%tev%XS_hjW_ratio(j)
!   XS_rat(2) = t%tev%XS_hjZ_ratio(j) 
      
!   BR_rat(1) = div( t%BR_hjWW(j) , t%BR_HWW_SM(j)  ,0.0D0,1.0D0) 
!   BR_rat(2) = div( t%BR_hjZZ(j) , t%BR_HZZ_SM(j)  ,0.0D0,1.0D0)

   channel_rat(1,:) = (/ t%tev%XS_hjW_ratio(j) , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(2,:) = (/ t%tev%XS_hjZ_ratio(j) , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(3,:) = (/ t%tev%XS_hjW_ratio(j) , div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(4,:) = (/ t%tev%XS_hjZ_ratio(j) , div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)

   channel_SM(1,:) = (/ t%tev%XS_HW_SM(j) , t%BR_HWW_SM(j) /) 
   channel_SM(2,:) = (/ t%tev%XS_HZ_SM(j) , t%BR_HWW_SM(j) /)
   channel_SM(3,:) = (/ t%tev%XS_HW_SM(j) , t%BR_HZZ_SM(j) /) 
   channel_SM(4,:) = (/ t%tev%XS_HZ_SM(j) , t%BR_HZZ_SM(j) /)
  

  case(6008,9998)
!   ns = 5; nb = 5; call initialise_XS_rat_BR_rat     
   nc = 25; call initialise_channel_rat_SM

!   XS_rat(1) = t%tev%XS_hjW_ratio(j)
!   XS_rat(2) = t%tev%XS_hj_ratio(j)  
!   XS_rat(3) = t%tev%XS_hjZ_ratio(j)  
!   XS_rat(4) = t%tev%XS_vbf_ratio(j)
!   XS_rat(5) = t%tev%XS_tthj_ratio(j)

!   BR_rat(1) = div(    t%BR_hjbb(j) , t%BR_Hbb_SM(j)     ,0.0D0,1.0D0)
!   BR_rat(2) = div(    t%BR_hjWW(j) , t%BR_HWW_SM(j)     ,0.0D0,1.0D0)   
!   BR_rat(3) = div(t%BR_hjtautau(j) , t%BR_Htautau_SM(j) ,0.0D0,1.0D0)
!   BR_rat(4) = div(t%BR_hjgaga(j) , t%BR_Hgaga_SM(j)  ,0.0D0,1.0D0)
!   BR_rat(5) = div(( t%BR_hjss(j)+t%BR_hjcc(j)+t%BR_hjbb(j)+t%BR_hjgg(j) ) &
!             & , t%BR_Hjets_SM(j) ,0.0D0,1.0D0)

   channel_rat(1,:) = (/ t%tev%XS_hj_ratio(j)  , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(2,:) = (/ t%tev%XS_vbf_ratio(j) , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(3,:) = (/ t%tev%XS_hjW_ratio(j) , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(4,:) = (/ t%tev%XS_hjZ_ratio(j) , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(5,:) = (/ t%tev%XS_tthj_ratio(j) , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(6,:) = (/ t%tev%XS_hj_ratio(j)  , div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(7,:) = (/ t%tev%XS_vbf_ratio(j) , div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(8,:) = (/ t%tev%XS_hjW_ratio(j) , div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(9,:) = (/ t%tev%XS_hjZ_ratio(j) , div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(10,:) = (/ t%tev%XS_tthj_ratio(j) , div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(11,:) = (/ t%tev%XS_hj_ratio(j)  , &
   & div((t%BR_hjss(j)+t%BR_hjcc(j)+t%BR_hjbb(j)+t%BR_hjgg(j)),t%BR_Hjets_SM(j),0.0D0,1.0D0) /)
   channel_rat(12,:) = (/ t%tev%XS_vbf_ratio(j) , &
   & div((t%BR_hjss(j)+t%BR_hjcc(j)+t%BR_hjbb(j)+t%BR_hjgg(j)),t%BR_Hjets_SM(j),0.0D0,1.0D0) /)
   channel_rat(13,:) = (/ t%tev%XS_hjW_ratio(j) , &
   & div((t%BR_hjss(j)+t%BR_hjcc(j)+t%BR_hjbb(j)+t%BR_hjgg(j)),t%BR_Hjets_SM(j),0.0D0,1.0D0) /)
   channel_rat(14,:) = (/ t%tev%XS_hjZ_ratio(j) , &
   & div((t%BR_hjss(j)+t%BR_hjcc(j)+t%BR_hjbb(j)+t%BR_hjgg(j)),t%BR_Hjets_SM(j),0.0D0,1.0D0) /)
   channel_rat(15,:) = (/ t%tev%XS_tthj_ratio(j) , &
   & div((t%BR_hjss(j)+t%BR_hjcc(j)+t%BR_hjbb(j)+t%BR_hjgg(j)),t%BR_Hjets_SM(j),0.0D0,1.0D0) /)
   channel_rat(16,:) = (/ t%tev%XS_hj_ratio(j)  , div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)
   channel_rat(17,:) = (/ t%tev%XS_vbf_ratio(j) , div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)
   channel_rat(18,:) = (/ t%tev%XS_hjW_ratio(j) , div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)
   channel_rat(19,:) = (/ t%tev%XS_hjZ_ratio(j) , div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)
   channel_rat(20,:) = (/ t%tev%XS_tthj_ratio(j) , div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)   
   channel_rat(21,:) = (/ t%tev%XS_hj_ratio(j)  , div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0) /)
   channel_rat(22,:) = (/ t%tev%XS_vbf_ratio(j) , div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0) /)
   channel_rat(23,:) = (/ t%tev%XS_hjW_ratio(j) , div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0) /)
   channel_rat(24,:) = (/ t%tev%XS_hjZ_ratio(j) , div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0) /)
   channel_rat(25,:) = (/ t%tev%XS_tthj_ratio(j) , div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0) /)
   
   channel_SM(1,:) = (/ t%tev%XS_H_SM(j)   , t%BR_HWW_SM(j) /)
   channel_SM(2,:) = (/ t%tev%XS_vbf_SM(j) , t%BR_HWW_SM(j) /)
   channel_SM(3,:) = (/ t%tev%XS_HW_SM(j) , t%BR_HWW_SM(j) /) 
   channel_SM(4,:) = (/ t%tev%XS_HZ_SM(j) , t%BR_HWW_SM(j) /)
   channel_SM(5,:) = (/ t%tev%XS_ttH_SM(j) , t%BR_HWW_SM(j) /)
   channel_SM(6,:) = (/ t%tev%XS_H_SM(j)   , t%BR_Htautau_SM(j) /)
   channel_SM(7,:) = (/ t%tev%XS_vbf_SM(j) , t%BR_Htautau_SM(j) /)
   channel_SM(8,:) = (/ t%tev%XS_HW_SM(j) , t%BR_Htautau_SM(j) /) 
   channel_SM(9,:) = (/ t%tev%XS_HZ_SM(j) , t%BR_Htautau_SM(j) /)
   channel_SM(10,:) = (/ t%tev%XS_ttH_SM(j) , t%BR_Htautau_SM(j) /)
   channel_SM(11,:) = (/ t%tev%XS_H_SM(j)   , t%BR_Hjets_SM(j) /)
   channel_SM(12,:) = (/ t%tev%XS_vbf_SM(j) , t%BR_Hjets_SM(j) /)
   channel_SM(13,:) = (/ t%tev%XS_HW_SM(j) , t%BR_Hjets_SM(j) /) 
   channel_SM(14,:) = (/ t%tev%XS_HZ_SM(j) , t%BR_Hjets_SM(j) /)
   channel_SM(15,:) = (/ t%tev%XS_ttH_SM(j) , t%BR_Hjets_SM(j) /)
   channel_SM(16,:) = (/ t%tev%XS_H_SM(j)   , t%BR_Hgaga_SM(j) /)
   channel_SM(17,:) = (/ t%tev%XS_vbf_SM(j) , t%BR_Hgaga_SM(j) /)
   channel_SM(18,:) = (/ t%tev%XS_HW_SM(j) , t%BR_Hgaga_SM(j) /) 
   channel_SM(19,:) = (/ t%tev%XS_HZ_SM(j) , t%BR_Hgaga_SM(j) /)
   channel_SM(20,:) = (/ t%tev%XS_ttH_SM(j) , t%BR_Hgaga_SM(j) /)
   channel_SM(21,:) = (/ t%tev%XS_H_SM(j)   , t%BR_Hbb_SM(j) /)
   channel_SM(22,:) = (/ t%tev%XS_vbf_SM(j) , t%BR_Hbb_SM(j) /)
   channel_SM(23,:) = (/ t%tev%XS_HW_SM(j) , t%BR_Hbb_SM(j) /) 
   channel_SM(24,:) = (/ t%tev%XS_HZ_SM(j) , t%BR_Hbb_SM(j) /)   
   channel_SM(25,:) = (/ t%tev%XS_ttH_SM(j) , t%BR_Hbb_SM(j) /)   

  case(6183,3233)
!   ns = 4; nb = 4; call initialise_XS_rat_BR_rat    
   nc = 16; call initialise_channel_rat_SM

!   XS_rat(1) = t%tev%XS_hjW_ratio(j)
!   XS_rat(2) = t%tev%XS_hj_ratio(j)  
!   XS_rat(3) = t%tev%XS_hjZ_ratio(j)  
!   XS_rat(4) = t%tev%XS_vbf_ratio(j)

!   BR_rat(1) = div(t%BR_hjWW(j) , t%BR_HWW_SM(j)         ,0.0D0,1.0D0)
!   BR_rat(2) = div( t%BR_hjZZ(j) , t%BR_HZZ_SM(j)        ,0.0D0,1.0D0)    
!   BR_rat(3) = div(t%BR_hjtautau(j) , t%BR_Htautau_SM(j) ,0.0D0,1.0D0)   
!   BR_rat(4) = div(t%BR_hjgaga(j) , t%BR_Hgaga_SM(j)     ,0.0D0,1.0D0)   

   channel_rat(1,:) = (/ t%tev%XS_hj_ratio(j)  , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(2,:) = (/ t%tev%XS_vbf_ratio(j) , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(3,:) = (/ t%tev%XS_hjW_ratio(j) , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(4,:) = (/ t%tev%XS_hjZ_ratio(j) , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(5,:) = (/ t%tev%XS_hj_ratio(j)  , div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(6,:) = (/ t%tev%XS_vbf_ratio(j) , div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(7,:) = (/ t%tev%XS_hjW_ratio(j) , div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(8,:) = (/ t%tev%XS_hjZ_ratio(j) , div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(9,:) = (/ t%tev%XS_hj_ratio(j)  , div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(10,:) = (/ t%tev%XS_vbf_ratio(j) , div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(11,:) = (/ t%tev%XS_hjW_ratio(j) , div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(12,:) = (/ t%tev%XS_hjZ_ratio(j) , div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(13,:) = (/ t%tev%XS_hj_ratio(j)  , div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)
   channel_rat(14,:) = (/ t%tev%XS_vbf_ratio(j) , div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)
   channel_rat(15,:) = (/ t%tev%XS_hjW_ratio(j) , div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)
   channel_rat(16,:) = (/ t%tev%XS_hjZ_ratio(j) , div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)

   channel_SM(1,:) = (/ t%tev%XS_H_SM(j)   , t%BR_HWW_SM(j) /)
   channel_SM(2,:) = (/ t%tev%XS_vbf_SM(j) , t%BR_HWW_SM(j) /)
   channel_SM(3,:) = (/ t%tev%XS_HW_SM(j) , t%BR_HWW_SM(j) /) 
   channel_SM(4,:) = (/ t%tev%XS_HZ_SM(j) , t%BR_HWW_SM(j) /)
   channel_SM(5,:) = (/ t%tev%XS_H_SM(j)   , t%BR_Htautau_SM(j) /)
   channel_SM(6,:) = (/ t%tev%XS_vbf_SM(j) , t%BR_Htautau_SM(j) /)
   channel_SM(7,:) = (/ t%tev%XS_HW_SM(j) , t%BR_Htautau_SM(j) /) 
   channel_SM(8,:) = (/ t%tev%XS_HZ_SM(j) , t%BR_Htautau_SM(j) /)
   channel_SM(9,:) = (/ t%tev%XS_H_SM(j)   , t%BR_HZZ_SM(j) /)
   channel_SM(10,:) = (/ t%tev%XS_vbf_SM(j) , t%BR_HZZ_SM(j) /)
   channel_SM(11,:) = (/ t%tev%XS_HW_SM(j) , t%BR_HZZ_SM(j) /) 
   channel_SM(12,:) = (/ t%tev%XS_HZ_SM(j) , t%BR_HZZ_SM(j) /)
   channel_SM(13,:) = (/ t%tev%XS_H_SM(j)   , t%BR_Hgaga_SM(j) /)
   channel_SM(14,:) = (/ t%tev%XS_vbf_SM(j) , t%BR_Hgaga_SM(j) /)
   channel_SM(15,:) = (/ t%tev%XS_HW_SM(j) , t%BR_Hgaga_SM(j) /) 
   channel_SM(16,:) = (/ t%tev%XS_HZ_SM(j) , t%BR_Hgaga_SM(j) /)

  case(6286)
!   ns = 4; nb = 4; call initialise_XS_rat_BR_rat    
   nc = 9; call initialise_channel_rat_SM
 
   channel_rat(1,:) = (/ t%tev%XS_hj_ratio(j)  , div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(2,:) = (/ t%tev%XS_vbf_ratio(j) , div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(3,:) = (/ t%tev%XS_hjW_ratio(j) , div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(4,:) = (/ t%tev%XS_hjW_ratio(j) , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(5,:) = (/ t%tev%XS_hjW_ratio(j) , div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(6,:) = (/ t%tev%XS_hjZ_ratio(j) , div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(7,:) = (/ t%tev%XS_hjZ_ratio(j) , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(8,:) = (/ t%tev%XS_hjZ_ratio(j) , div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(9,:) = (/ t%tev%XS_hjZ_ratio(j) , div(t%BR_hjmumu(j),t%BR_Hmumu_SM(j),0.0D0,1.0D0) /)

   channel_SM(1,:) = (/ t%tev%XS_H_SM(j)  , t%BR_HZZ_SM(j) /)
   channel_SM(2,:) = (/ t%tev%XS_vbf_SM(j), t%BR_HZZ_SM(j) /)
   channel_SM(3,:) = (/ t%tev%XS_HW_SM(j) , t%BR_HZZ_SM(j) /) 
   channel_SM(4,:) = (/ t%tev%XS_HW_SM(j) , t%BR_HWW_SM(j) /) 
   channel_SM(5,:) = (/ t%tev%XS_HW_SM(j) , t%BR_Htautau_SM(j) /) 
   channel_SM(6,:) = (/ t%tev%XS_HZ_SM(j) , t%BR_HZZ_SM(j) /)
   channel_SM(7,:) = (/ t%tev%XS_HZ_SM(j) , t%BR_HWW_SM(j) /)
   channel_SM(8,:) = (/ t%tev%XS_HZ_SM(j) , t%BR_Htautau_SM(j) /)
   channel_SM(9,:) = (/ t%tev%XS_HZ_SM(j) , t%BR_Hmumu_SM(j) /)

  case(6305)
!   ns = 4; nb = 3; call initialise_XS_rat_BR_rat    
   nc = 12; call initialise_channel_rat_SM

!   XS_rat(1) = t%tev%XS_hjW_ratio(j)
!   XS_rat(2) = t%tev%XS_hj_ratio(j)  
!   XS_rat(3) = t%tev%XS_hjZ_ratio(j)  
!   XS_rat(4) = t%tev%XS_vbf_ratio(j)

!   BR_rat(1) = div(t%BR_hjWW(j) , t%BR_HWW_SM(j)         ,0.0D0,1.0D0)
!   BR_rat(2) = div( t%BR_hjZZ(j) , t%BR_HZZ_SM(j)        ,0.0D0,1.0D0)    
!   BR_rat(3) = div(t%BR_hjtautau(j) , t%BR_Htautau_SM(j) ,0.0D0,1.0D0)     

   channel_rat(1,:) = (/ t%tev%XS_hj_ratio(j)  , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(2,:) = (/ t%tev%XS_vbf_ratio(j) , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(3,:) = (/ t%tev%XS_hjW_ratio(j) , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(4,:) = (/ t%tev%XS_hjZ_ratio(j) , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(5,:) = (/ t%tev%XS_hj_ratio(j)  , div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(6,:) = (/ t%tev%XS_vbf_ratio(j) , div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(7,:) = (/ t%tev%XS_hjW_ratio(j) , div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(8,:) = (/ t%tev%XS_hjZ_ratio(j) , div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(9,:) = (/ t%tev%XS_hj_ratio(j)  , div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(10,:) = (/ t%tev%XS_vbf_ratio(j) , div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(11,:) = (/ t%tev%XS_hjW_ratio(j) , div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(12,:) = (/ t%tev%XS_hjZ_ratio(j) , div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)

   channel_SM(1,:) = (/ t%tev%XS_H_SM(j)   , t%BR_HWW_SM(j) /)
   channel_SM(2,:) = (/ t%tev%XS_vbf_SM(j) , t%BR_HWW_SM(j) /)
   channel_SM(3,:) = (/ t%tev%XS_HW_SM(j) , t%BR_HWW_SM(j) /) 
   channel_SM(4,:) = (/ t%tev%XS_HZ_SM(j) , t%BR_HWW_SM(j) /)
   channel_SM(5,:) = (/ t%tev%XS_H_SM(j)   , t%BR_Htautau_SM(j) /)
   channel_SM(6,:) = (/ t%tev%XS_vbf_SM(j) , t%BR_Htautau_SM(j) /)
   channel_SM(7,:) = (/ t%tev%XS_HW_SM(j) , t%BR_Htautau_SM(j) /) 
   channel_SM(8,:) = (/ t%tev%XS_HZ_SM(j) , t%BR_Htautau_SM(j) /)
   channel_SM(9,:) = (/ t%tev%XS_H_SM(j)   , t%BR_HZZ_SM(j) /)
   channel_SM(10,:) = (/ t%tev%XS_vbf_SM(j) , t%BR_HZZ_SM(j) /)
   channel_SM(11,:) = (/ t%tev%XS_HW_SM(j) , t%BR_HZZ_SM(j) /) 
   channel_SM(12,:) = (/ t%tev%XS_HZ_SM(j) , t%BR_HZZ_SM(j) /)

  case(6096,10606,10806,10884)
   nc = 30; call initialise_channel_rat_SM

   channel_rat(1,:) = (/ t%tev%XS_hj_ratio(j)  , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(2,:) = (/ t%tev%XS_vbf_ratio(j) , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(3,:) = (/ t%tev%XS_hjW_ratio(j) , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(4,:) = (/ t%tev%XS_hjZ_ratio(j) , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(5,:) = (/ t%tev%XS_tthj_ratio(j) , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(6,:) = (/ t%tev%XS_hj_ratio(j)  , div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(7,:) = (/ t%tev%XS_vbf_ratio(j) , div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(8,:) = (/ t%tev%XS_hjW_ratio(j) , div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(9,:) = (/ t%tev%XS_hjZ_ratio(j) , div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(10,:) = (/ t%tev%XS_tthj_ratio(j) , div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(11,:) = (/ t%tev%XS_hj_ratio(j)  , &
   & div((t%BR_hjss(j)+t%BR_hjcc(j)+t%BR_hjbb(j)+t%BR_hjgg(j)),t%BR_Hjets_SM(j),0.0D0,1.0D0) /)
   channel_rat(12,:) = (/ t%tev%XS_vbf_ratio(j) , &
   & div((t%BR_hjss(j)+t%BR_hjcc(j)+t%BR_hjbb(j)+t%BR_hjgg(j)),t%BR_Hjets_SM(j),0.0D0,1.0D0) /)
   channel_rat(13,:) = (/ t%tev%XS_hjW_ratio(j) , &
   & div((t%BR_hjss(j)+t%BR_hjcc(j)+t%BR_hjbb(j)+t%BR_hjgg(j)),t%BR_Hjets_SM(j),0.0D0,1.0D0) /)
   channel_rat(14,:) = (/ t%tev%XS_hjZ_ratio(j) , &
   & div((t%BR_hjss(j)+t%BR_hjcc(j)+t%BR_hjbb(j)+t%BR_hjgg(j)),t%BR_Hjets_SM(j),0.0D0,1.0D0) /)
   channel_rat(15,:) = (/ t%tev%XS_tthj_ratio(j) , &
   & div((t%BR_hjss(j)+t%BR_hjcc(j)+t%BR_hjbb(j)+t%BR_hjgg(j)),t%BR_Hjets_SM(j),0.0D0,1.0D0) /)
   channel_rat(16,:) = (/ t%tev%XS_hj_ratio(j)  , div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(17,:) = (/ t%tev%XS_vbf_ratio(j) , div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(18,:) = (/ t%tev%XS_hjW_ratio(j) , div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(19,:) = (/ t%tev%XS_hjZ_ratio(j) , div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(20,:) = (/ t%tev%XS_tthj_ratio(j) , div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)      
   channel_rat(21,:) = (/ t%tev%XS_hj_ratio(j)  , div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)
   channel_rat(22,:) = (/ t%tev%XS_vbf_ratio(j) , div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)
   channel_rat(23,:) = (/ t%tev%XS_hjW_ratio(j) , div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)
   channel_rat(24,:) = (/ t%tev%XS_hjZ_ratio(j) , div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)
   channel_rat(25,:) = (/ t%tev%XS_tthj_ratio(j) , div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)   
   channel_rat(26,:) = (/ t%tev%XS_hj_ratio(j)  , div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0) /)
   channel_rat(27,:) = (/ t%tev%XS_vbf_ratio(j) , div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0) /)
   channel_rat(28,:) = (/ t%tev%XS_hjW_ratio(j) , div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0) /)
   channel_rat(29,:) = (/ t%tev%XS_hjZ_ratio(j) , div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0) /)
   channel_rat(30,:) = (/ t%tev%XS_tthj_ratio(j) , div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0) /)
   
   channel_SM(1,:) = (/ t%tev%XS_H_SM(j)   , t%BR_HWW_SM(j) /)
   channel_SM(2,:) = (/ t%tev%XS_vbf_SM(j) , t%BR_HWW_SM(j) /)
   channel_SM(3,:) = (/ t%tev%XS_HW_SM(j) , t%BR_HWW_SM(j) /) 
   channel_SM(4,:) = (/ t%tev%XS_HZ_SM(j) , t%BR_HWW_SM(j) /)
   channel_SM(5,:) = (/ t%tev%XS_ttH_SM(j) , t%BR_HWW_SM(j) /)
   channel_SM(6,:) = (/ t%tev%XS_H_SM(j)   , t%BR_Htautau_SM(j) /)
   channel_SM(7,:) = (/ t%tev%XS_vbf_SM(j) , t%BR_Htautau_SM(j) /)
   channel_SM(8,:) = (/ t%tev%XS_HW_SM(j) , t%BR_Htautau_SM(j) /) 
   channel_SM(9,:) = (/ t%tev%XS_HZ_SM(j) , t%BR_Htautau_SM(j) /)
   channel_SM(10,:) = (/ t%tev%XS_ttH_SM(j) , t%BR_Htautau_SM(j) /)
   channel_SM(11,:) = (/ t%tev%XS_H_SM(j)   , t%BR_Hjets_SM(j) /)
   channel_SM(12,:) = (/ t%tev%XS_vbf_SM(j) , t%BR_Hjets_SM(j) /)
   channel_SM(13,:) = (/ t%tev%XS_HW_SM(j) , t%BR_Hjets_SM(j) /) 
   channel_SM(14,:) = (/ t%tev%XS_HZ_SM(j) , t%BR_Hjets_SM(j) /)
   channel_SM(15,:) = (/ t%tev%XS_ttH_SM(j) , t%BR_Hjets_SM(j) /)
   channel_SM(16,:) = (/ t%tev%XS_H_SM(j)   , t%BR_HZZ_SM(j) /)
   channel_SM(17,:) = (/ t%tev%XS_vbf_SM(j) , t%BR_HZZ_SM(j) /)
   channel_SM(18,:) = (/ t%tev%XS_HW_SM(j) , t%BR_HZZ_SM(j) /) 
   channel_SM(19,:) = (/ t%tev%XS_HZ_SM(j) , t%BR_HZZ_SM(j) /)
   channel_SM(20,:) = (/ t%tev%XS_ttH_SM(j) , t%BR_HZZ_SM(j) /)
   channel_SM(21,:) = (/ t%tev%XS_H_SM(j)   , t%BR_Hgaga_SM(j) /)
   channel_SM(22,:) = (/ t%tev%XS_vbf_SM(j) , t%BR_Hgaga_SM(j) /)
   channel_SM(23,:) = (/ t%tev%XS_HW_SM(j) , t%BR_Hgaga_SM(j) /) 
   channel_SM(24,:) = (/ t%tev%XS_HZ_SM(j) , t%BR_Hgaga_SM(j) /)
   channel_SM(25,:) = (/ t%tev%XS_ttH_SM(j) , t%BR_Hgaga_SM(j) /)
   channel_SM(26,:) = (/ t%tev%XS_H_SM(j)   , t%BR_Hbb_SM(j) /)
   channel_SM(27,:) = (/ t%tev%XS_vbf_SM(j) , t%BR_Hbb_SM(j) /)
   channel_SM(28,:) = (/ t%tev%XS_HW_SM(j) , t%BR_Hbb_SM(j) /) 
   channel_SM(29,:) = (/ t%tev%XS_HZ_SM(j) , t%BR_Hbb_SM(j) /)   
   channel_SM(30,:) = (/ t%tev%XS_ttH_SM(j) , t%BR_Hbb_SM(j) /)   

  case(6229)
!   ns = 4; nb = 6; call initialise_XS_rat_BR_rat
   nc = 24; call initialise_channel_rat_SM

!   XS_rat(1) = t%tev%XS_hjW_ratio(j)
!   XS_rat(2) = t%tev%XS_hj_ratio(j)  
!   XS_rat(3) = t%tev%XS_hjZ_ratio(j)  
!   XS_rat(4) = t%tev%XS_vbf_ratio(j)

!   BR_rat(1) = div(    t%BR_hjbb(j) , t%BR_Hbb_SM(j)     ,0.0D0,1.0D0) 
!   BR_rat(2) = div(    t%BR_hjWW(j) , t%BR_HWW_SM(j)     ,0.0D0,1.0D0)   
!   BR_rat(3) = div(t%BR_hjtautau(j) , t%BR_Htautau_SM(j) ,0.0D0,1.0D0) 
!   BR_rat(4) = div(  t%BR_hjgaga(j) , t%BR_Hgaga_SM(j)   ,0.0D0,1.0D0) 
!   BR_rat(5) = div(( t%BR_hjss(j)+t%BR_hjcc(j)+t%BR_hjbb(j)+t%BR_hjgg(j) ) &
!             & , t%BR_Hjets_SM(j) ,0.0D0,1.0D0) 
!   BR_rat(6) = div(    t%BR_hjZZ(j) , t%BR_HZZ_SM(j)     ,0.0D0,1.0D0) 

   channel_rat(1,:) = (/ t%tev%XS_hj_ratio(j)  , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(2,:) = (/ t%tev%XS_vbf_ratio(j) , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(3,:) = (/ t%tev%XS_hjW_ratio(j) , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(4,:) = (/ t%tev%XS_hjZ_ratio(j) , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(5,:) = (/ t%tev%XS_hj_ratio(j)  , div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(6,:) = (/ t%tev%XS_vbf_ratio(j) , div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(7,:) = (/ t%tev%XS_hjW_ratio(j) , div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(8,:) = (/ t%tev%XS_hjZ_ratio(j) , div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(9,:) = (/ t%tev%XS_hj_ratio(j)  , &
   & div((t%BR_hjss(j)+t%BR_hjcc(j)+t%BR_hjbb(j)+t%BR_hjgg(j)),t%BR_Hjets_SM(j),0.0D0,1.0D0) /)
   channel_rat(10,:) = (/ t%tev%XS_vbf_ratio(j) , &
   & div((t%BR_hjss(j)+t%BR_hjcc(j)+t%BR_hjbb(j)+t%BR_hjgg(j)),t%BR_Hjets_SM(j),0.0D0,1.0D0) /)
   channel_rat(11,:) = (/ t%tev%XS_hjW_ratio(j) , &
   & div((t%BR_hjss(j)+t%BR_hjcc(j)+t%BR_hjbb(j)+t%BR_hjgg(j)),t%BR_Hjets_SM(j),0.0D0,1.0D0) /)
   channel_rat(12,:) = (/ t%tev%XS_hjZ_ratio(j) , &
   & div((t%BR_hjss(j)+t%BR_hjcc(j)+t%BR_hjbb(j)+t%BR_hjgg(j)),t%BR_Hjets_SM(j),0.0D0,1.0D0) /)
   channel_rat(13,:) = (/ t%tev%XS_hj_ratio(j)  , div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(14,:) = (/ t%tev%XS_vbf_ratio(j) , div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(15,:) = (/ t%tev%XS_hjW_ratio(j) , div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(16,:) = (/ t%tev%XS_hjZ_ratio(j) , div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(17,:) = (/ t%tev%XS_hj_ratio(j)  , div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)
   channel_rat(18,:) = (/ t%tev%XS_vbf_ratio(j) , div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)
   channel_rat(19,:) = (/ t%tev%XS_hjW_ratio(j) , div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)
   channel_rat(20,:) = (/ t%tev%XS_hjZ_ratio(j) , div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)
   channel_rat(21,:) = (/ t%tev%XS_hj_ratio(j)  , div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0) /)
   channel_rat(22,:) = (/ t%tev%XS_vbf_ratio(j) , div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0) /)
   channel_rat(23,:) = (/ t%tev%XS_hjW_ratio(j) , div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0) /)
   channel_rat(24,:) = (/ t%tev%XS_hjZ_ratio(j) , div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0) /)

   channel_SM(1,:) = (/ t%tev%XS_H_SM(j)   , t%BR_HWW_SM(j) /)
   channel_SM(2,:) = (/ t%tev%XS_vbf_SM(j) , t%BR_HWW_SM(j) /)
   channel_SM(3,:) = (/ t%tev%XS_HW_SM(j) , t%BR_HWW_SM(j) /) 
   channel_SM(4,:) = (/ t%tev%XS_HZ_SM(j) , t%BR_HWW_SM(j) /)
   channel_SM(5,:) = (/ t%tev%XS_H_SM(j)   , t%BR_Htautau_SM(j) /)
   channel_SM(6,:) = (/ t%tev%XS_vbf_SM(j) , t%BR_Htautau_SM(j) /)
   channel_SM(7,:) = (/ t%tev%XS_HW_SM(j) , t%BR_Htautau_SM(j) /) 
   channel_SM(8,:) = (/ t%tev%XS_HZ_SM(j) , t%BR_Htautau_SM(j) /)
   channel_SM(9,:) = (/ t%tev%XS_H_SM(j)   , t%BR_Hjets_SM(j) /)
   channel_SM(10,:) = (/ t%tev%XS_vbf_SM(j) , t%BR_Hjets_SM(j) /)
   channel_SM(11,:) = (/ t%tev%XS_HW_SM(j) , t%BR_Hjets_SM(j) /) 
   channel_SM(12,:) = (/ t%tev%XS_HZ_SM(j) , t%BR_Hjets_SM(j) /)
   channel_SM(13,:) = (/ t%tev%XS_H_SM(j)   , t%BR_HZZ_SM(j) /)
   channel_SM(14,:) = (/ t%tev%XS_vbf_SM(j) , t%BR_HZZ_SM(j) /)
   channel_SM(15,:) = (/ t%tev%XS_HW_SM(j) , t%BR_HZZ_SM(j) /) 
   channel_SM(16,:) = (/ t%tev%XS_HZ_SM(j) , t%BR_HZZ_SM(j) /)
   channel_SM(17,:) = (/ t%tev%XS_H_SM(j)   , t%BR_Hgaga_SM(j) /)
   channel_SM(18,:) = (/ t%tev%XS_vbf_SM(j) , t%BR_Hgaga_SM(j) /)
   channel_SM(19,:) = (/ t%tev%XS_HW_SM(j) , t%BR_Hgaga_SM(j) /) 
   channel_SM(20,:) = (/ t%tev%XS_HZ_SM(j) , t%BR_Hgaga_SM(j) /)
   channel_SM(21,:) = (/ t%tev%XS_H_SM(j)   , t%BR_Hbb_SM(j) /)
   channel_SM(22,:) = (/ t%tev%XS_vbf_SM(j) , t%BR_Hbb_SM(j) /)
   channel_SM(23,:) = (/ t%tev%XS_HW_SM(j) , t%BR_Hbb_SM(j) /) 
   channel_SM(24,:) = (/ t%tev%XS_HZ_SM(j) , t%BR_Hbb_SM(j) /)   
   
  case(6304)
!   ns = 4; nb = 5; call initialise_XS_rat_BR_rat
   nc = 20; call initialise_channel_rat_SM

   channel_rat(1,:) = (/ t%tev%XS_hj_ratio(j)  , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(2,:) = (/ t%tev%XS_vbf_ratio(j) , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(3,:) = (/ t%tev%XS_hjW_ratio(j) , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(4,:) = (/ t%tev%XS_hjZ_ratio(j) , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(5,:) = (/ t%tev%XS_hj_ratio(j)  , div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(6,:) = (/ t%tev%XS_vbf_ratio(j) , div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(7,:) = (/ t%tev%XS_hjW_ratio(j) , div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(8,:) = (/ t%tev%XS_hjZ_ratio(j) , div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(9,:) = (/ t%tev%XS_hj_ratio(j)  , div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(10,:) = (/ t%tev%XS_vbf_ratio(j) , div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(11,:) = (/ t%tev%XS_hjW_ratio(j) , div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(12,:) = (/ t%tev%XS_hjZ_ratio(j) , div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(13,:) = (/ t%tev%XS_hj_ratio(j)  , div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)
   channel_rat(14,:) = (/ t%tev%XS_vbf_ratio(j) , div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)
   channel_rat(15,:) = (/ t%tev%XS_hjW_ratio(j) , div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)
   channel_rat(16,:) = (/ t%tev%XS_hjZ_ratio(j) , div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)
   channel_rat(17,:) = (/ t%tev%XS_hj_ratio(j)  , div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0) /)
   channel_rat(18,:) = (/ t%tev%XS_vbf_ratio(j) , div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0) /)
   channel_rat(19,:) = (/ t%tev%XS_hjW_ratio(j) , div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0) /)
   channel_rat(20,:) = (/ t%tev%XS_hjZ_ratio(j) , div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0) /)

   channel_SM(1,:) = (/ t%tev%XS_H_SM(j)   , t%BR_HWW_SM(j) /)
   channel_SM(2,:) = (/ t%tev%XS_vbf_SM(j) , t%BR_HWW_SM(j) /)
   channel_SM(3,:) = (/ t%tev%XS_HW_SM(j) , t%BR_HWW_SM(j) /) 
   channel_SM(4,:) = (/ t%tev%XS_HZ_SM(j) , t%BR_HWW_SM(j) /)
   channel_SM(5,:) = (/ t%tev%XS_H_SM(j)   , t%BR_Htautau_SM(j) /)
   channel_SM(6,:) = (/ t%tev%XS_vbf_SM(j) , t%BR_Htautau_SM(j) /)
   channel_SM(7,:) = (/ t%tev%XS_HW_SM(j) , t%BR_Htautau_SM(j) /) 
   channel_SM(8,:) = (/ t%tev%XS_HZ_SM(j) , t%BR_Htautau_SM(j) /)
   channel_SM(9,:) = (/ t%tev%XS_H_SM(j)   , t%BR_HZZ_SM(j) /)
   channel_SM(10,:) = (/ t%tev%XS_vbf_SM(j) , t%BR_HZZ_SM(j) /)
   channel_SM(11,:) = (/ t%tev%XS_HW_SM(j) , t%BR_HZZ_SM(j) /) 
   channel_SM(12,:) = (/ t%tev%XS_HZ_SM(j) , t%BR_HZZ_SM(j) /)
   channel_SM(13,:) = (/ t%tev%XS_H_SM(j)   , t%BR_Hgaga_SM(j) /)
   channel_SM(14,:) = (/ t%tev%XS_vbf_SM(j) , t%BR_Hgaga_SM(j) /)
   channel_SM(15,:) = (/ t%tev%XS_HW_SM(j) , t%BR_Hgaga_SM(j) /) 
   channel_SM(16,:) = (/ t%tev%XS_HZ_SM(j) , t%BR_Hgaga_SM(j) /)
   channel_SM(17,:) = (/ t%tev%XS_H_SM(j)   , t%BR_Hbb_SM(j) /)
   channel_SM(18,:) = (/ t%tev%XS_vbf_SM(j) , t%BR_Hbb_SM(j) /)
   channel_SM(19,:) = (/ t%tev%XS_HW_SM(j) , t%BR_Hbb_SM(j) /) 
   channel_SM(20,:) = (/ t%tev%XS_HZ_SM(j) , t%BR_Hbb_SM(j) /)   
   
  case(6171)
!   ns = 4; nb = 3; call initialise_XS_rat_BR_rat
   nc = 12; call initialise_channel_rat_SM

!   XS_rat(1) = t%tev%XS_hjW_ratio(j)
!   XS_rat(2) = t%tev%XS_hj_ratio(j)  
!   XS_rat(3) = t%tev%XS_hjZ_ratio(j)  
!   XS_rat(4) = t%tev%XS_vbf_ratio(j)

!   BR_rat(1) = div(    t%BR_hjWW(j) , t%BR_HWW_SM(j)    ,0.0D0,1.0D0)   
!   BR_rat(2) = div(t%BR_hjtautau(j) , t%BR_Htautau_SM(j),0.0D0,1.0D0)  
!   BR_rat(3) = div(( t%BR_hjss(j)+t%BR_hjcc(j)+t%BR_hjbb(j)+t%BR_hjgg(j) ) &
!             & , t%BR_Hjets_SM(j),0.0D0,1.0D0) 

   channel_rat(1,:) = (/ t%tev%XS_hj_ratio(j)  , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(2,:) = (/ t%tev%XS_vbf_ratio(j) , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(3,:) = (/ t%tev%XS_hjW_ratio(j) , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(4,:) = (/ t%tev%XS_hjZ_ratio(j) , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(5,:) = (/ t%tev%XS_hj_ratio(j)  , div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(6,:) = (/ t%tev%XS_vbf_ratio(j) , div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(7,:) = (/ t%tev%XS_hjW_ratio(j) , div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(8,:) = (/ t%tev%XS_hjZ_ratio(j) , div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(9,:) = (/ t%tev%XS_hj_ratio(j)  , &
   & div((t%BR_hjss(j)+t%BR_hjcc(j)+t%BR_hjbb(j)+t%BR_hjgg(j)),t%BR_Hjets_SM(j),0.0D0,1.0D0) /)
   channel_rat(10,:) = (/ t%tev%XS_vbf_ratio(j) , &
   & div((t%BR_hjss(j)+t%BR_hjcc(j)+t%BR_hjbb(j)+t%BR_hjgg(j)),t%BR_Hjets_SM(j),0.0D0,1.0D0) /)
   channel_rat(11,:) = (/ t%tev%XS_hjW_ratio(j) , &
   & div((t%BR_hjss(j)+t%BR_hjcc(j)+t%BR_hjbb(j)+t%BR_hjgg(j)),t%BR_Hjets_SM(j),0.0D0,1.0D0) /)
   channel_rat(12,:) = (/ t%tev%XS_hjZ_ratio(j) , &
   & div((t%BR_hjss(j)+t%BR_hjcc(j)+t%BR_hjbb(j)+t%BR_hjgg(j)),t%BR_Hjets_SM(j),0.0D0,1.0D0) /)

   channel_SM(1,:) = (/ t%tev%XS_H_SM(j)   , t%BR_HWW_SM(j) /)
   channel_SM(2,:) = (/ t%tev%XS_vbf_SM(j) , t%BR_HWW_SM(j) /)
   channel_SM(3,:) = (/ t%tev%XS_HW_SM(j) , t%BR_HWW_SM(j) /) 
   channel_SM(4,:) = (/ t%tev%XS_HZ_SM(j) , t%BR_HWW_SM(j) /)
   channel_SM(5,:) = (/ t%tev%XS_H_SM(j)   , t%BR_Htautau_SM(j) /)
   channel_SM(6,:) = (/ t%tev%XS_vbf_SM(j) , t%BR_Htautau_SM(j) /)
   channel_SM(7,:) = (/ t%tev%XS_HW_SM(j) , t%BR_Htautau_SM(j) /) 
   channel_SM(8,:) = (/ t%tev%XS_HZ_SM(j) , t%BR_Htautau_SM(j) /)
   channel_SM(9,:) = (/ t%tev%XS_H_SM(j)   , t%BR_Hjets_SM(j) /)
   channel_SM(10,:) = (/ t%tev%XS_vbf_SM(j) , t%BR_Hjets_SM(j) /)
   channel_SM(11,:) = (/ t%tev%XS_HW_SM(j) , t%BR_Hjets_SM(j) /) 
   channel_SM(12,:) = (/ t%tev%XS_HZ_SM(j) , t%BR_Hjets_SM(j) /)

!---------------------- LHC 7/8 TeV searches --------------------  
  case(5757,2011005,3615,2012018)
!   ns = 2; nb = 1; call initialise_XS_rat_BR_rat 
   nc = 2; call initialise_channel_rat_SM
        
!   XS_rat(1) = t%tev%XS_hj_ratio(j) 
!   XS_rat(2) = t%tev%XS_vbf_ratio(j)
  
!   BR_rat(1) = div( t%BR_hjWW(j) , t%BR_HWW_SM(j)  ,0.0D0,1.0D0)

   channel_rat(1,:) = (/ t%lhc7%XS_hj_ratio(j) , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(2,:) = (/ t%lhc7%XS_vbf_ratio(j), div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)

   channel_SM(1,:) = (/ t%lhc7%XS_H_SM(j) , t%BR_HWW_SM(j) /)
   channel_SM(2,:) = (/ t%lhc7%XS_vbf_SM(j), t%BR_HWW_SM(j) /)

  case(12046)
   nc = 2; call initialise_channel_rat_SM
        
   channel_rat(1,:) = (/ t%lhc8%XS_hj_ratio(j) , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(2,:) = (/ t%lhc8%XS_vbf_ratio(j), div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)

   channel_SM(1,:) = (/ t%lhc8%XS_H_SM(j) , t%BR_HWW_SM(j) /)
   channel_SM(2,:) = (/ t%lhc8%XS_vbf_SM(j), t%BR_HWW_SM(j) /)


  case(2011048,11004,11015,11013,11028,11006,11017,110271,110272,14161,14162,5064,2012017,2011150,2011131) 
!   ns = 2; nb = 1; call initialise_XS_rat_BR_rat         
   nc = 2; call initialise_channel_rat_SM
   
!   XS_rat(1) = t%lhc7%XS_hj_ratio(j)   
!   XS_rat(2) = t%lhc7%XS_vbf_ratio(j)  
!   BR_rat(1) = div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0)

   channel_rat(1,:) = (/ t%lhc7%XS_hj_ratio(j) , div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(2,:) = (/ t%lhc7%XS_vbf_ratio(j) , div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /) 
   
   channel_SM(1,:) = (/ t%lhc7%XS_H_SM(j) , t%BR_HZZ_SM(j) /)
   channel_SM(2,:) = (/ t%lhc7%XS_vbf_SM(j) , t%BR_HZZ_SM(j) /)      

  case(11034,12039,13009,12006,2012078) 
!   ns = 1; nb = 2; call initialise_XS_rat_BR_rat         
   nc = 2; call initialise_channel_rat_SM
   
   channel_rat(1,:) = (/ t%lhc7%XS_hjW_ratio(j) , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(2,:) = (/ t%lhc7%XS_hjW_ratio(j) , div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /) 
   
   channel_SM(1,:) = (/ t%lhc7%XS_HW_SM(j) , t%BR_HWW_SM(j) /)
   channel_SM(2,:) = (/ t%lhc7%XS_HW_SM(j) , t%BR_Htautau_SM(j) /)  

  case(12051) 
   nc = 2; call initialise_channel_rat_SM
   
   channel_rat(1,:) = (/ t%lhc7%XS_hjW_ratio(j) , div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(2,:) = (/ t%lhc7%XS_hjZ_ratio(j) , div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /) 
   
   channel_SM(1,:) = (/ t%lhc7%XS_HW_SM(j) , t%BR_Htautau_SM(j) /)  
   channel_SM(2,:) = (/ t%lhc7%XS_HZ_SM(j) , t%BR_Htautau_SM(j) /)  


  case(2012015)
!   ns = 2; nb = 1; call initialise_XS_rat_BR_rat      
   nc = 2; call initialise_channel_rat_SM
  
!   XS_rat(1) = t%tev%XS_hjW_ratio(j)
!   XS_rat(2) = t%tev%XS_hjZ_ratio(j)  
  
!   BR_rat(1) = div(t%BR_hjbb(j) , t%BR_Hbb_SM(j),0.0D0,1.0D0)
 
   channel_rat(1,:) = (/ t%lhc7%XS_hjW_ratio(j) , div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0) /)
   channel_rat(2,:) = (/ t%lhc7%XS_hjZ_ratio(j) , div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0) /) 

   channel_SM(1,:) = (/ t%lhc7%XS_HW_SM(j) , t%BR_Hbb_SM(j) /) 
   channel_SM(2,:) = (/ t%lhc7%XS_HZ_SM(j) , t%BR_Hbb_SM(j) /)  

  case(2011162,1415) 
!   ns = 4; nb = 1; call initialise_XS_rat_BR_rat   
   nc = 4; call initialise_channel_rat_SM
        
!   XS_rat(1) = t%lhc7%XS_hj_ratio(j)   
!   XS_rat(2) = t%lhc7%XS_vbf_ratio(j)
!   XS_rat(3) = t%lhc7%XS_hjZ_ratio(j)  
!   XS_rat(4) = t%lhc7%XS_hjW_ratio(j) 
  
!   BR_rat(1) = div(t%BR_hjZZ(j)  , t%BR_HZZ_SM(j) ,0.0D0,1.0D0)

   channel_rat(1,:) = (/ t%lhc7%XS_hj_ratio(j) , div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(2,:) = (/ t%lhc7%XS_vbf_ratio(j), div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(3,:) = (/ t%lhc7%XS_hjZ_ratio(j), div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(4,:) = (/ t%lhc7%XS_hjW_ratio(j), div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)

   channel_SM(1,:) = (/ t%lhc7%XS_H_SM(j) , t%BR_HZZ_SM(j) /)
   channel_SM(2,:) = (/ t%lhc7%XS_vbf_SM(j), t%BR_HZZ_SM(j) /)
   channel_SM(3,:) = (/ t%lhc7%XS_HZ_SM(j), t%BR_HZZ_SM(j) /)
   channel_SM(4,:) = (/ t%lhc7%XS_HW_SM(j), t%BR_HZZ_SM(j) /)

  case(2012092,20130131) 
   nc = 4; call initialise_channel_rat_SM
        
   channel_rat(1,:) = (/ t%lhc8%XS_hj_ratio(j) , div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(2,:) = (/ t%lhc8%XS_vbf_ratio(j), div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(3,:) = (/ t%lhc8%XS_hjZ_ratio(j), div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(4,:) = (/ t%lhc8%XS_hjW_ratio(j), div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)

   channel_SM(1,:) = (/ t%lhc8%XS_H_SM(j) , t%BR_HZZ_SM(j) /)
   channel_SM(2,:) = (/ t%lhc8%XS_vbf_SM(j), t%BR_HZZ_SM(j) /)
   channel_SM(3,:) = (/ t%lhc8%XS_HZ_SM(j), t%BR_HZZ_SM(j) /)
   channel_SM(4,:) = (/ t%lhc8%XS_HW_SM(j), t%BR_HZZ_SM(j) /)

  case(11025,1997) 
!   ns = 5; nb = 1; call initialise_XS_rat_BR_rat   
   nc = 5; call initialise_channel_rat_SM
        
!   XS_rat(1) = t%lhc7%XS_hj_ratio(j)   
!   XS_rat(2) = t%lhc7%XS_vbf_ratio(j)
!   XS_rat(3) = t%lhc7%XS_hjZ_ratio(j)  
!   XS_rat(4) = t%lhc7%XS_hjW_ratio(j) 
!   XS_rat(5) = t%lhc7%XS_tthj_ratio(j) 
     
!   BR_rat(1) = div(t%BR_hjZZ(j)  , t%BR_HZZ_SM(j) ,0.0D0,1.0D0)

   channel_rat(1,:) = (/ t%lhc7%XS_hj_ratio(j) , div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(2,:) = (/ t%lhc7%XS_vbf_ratio(j), div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(3,:) = (/ t%lhc7%XS_hjZ_ratio(j), div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(4,:) = (/ t%lhc7%XS_hjW_ratio(j), div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(5,:) = (/ t%lhc7%XS_tthj_ratio(j), div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)

   channel_SM(1,:) = (/ t%lhc7%XS_H_SM(j) , t%BR_HZZ_SM(j) /)
   channel_SM(2,:) = (/ t%lhc7%XS_vbf_SM(j), t%BR_HZZ_SM(j) /)
   channel_SM(3,:) = (/ t%lhc7%XS_HZ_SM(j), t%BR_HZZ_SM(j) /)
   channel_SM(4,:) = (/ t%lhc7%XS_HW_SM(j), t%BR_HZZ_SM(j) /)
   channel_SM(5,:) = (/ t%lhc7%XS_ttH_SM(j), t%BR_HZZ_SM(j) /)

  case(12041,130021,130022)
   nc = 5; call initialise_channel_rat_SM
 
   channel_rat(1,:) = (/ t%lhc8%XS_hj_ratio(j) , div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(2,:) = (/ t%lhc8%XS_vbf_ratio(j), div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(3,:) = (/ t%lhc8%XS_hjZ_ratio(j), div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(4,:) = (/ t%lhc8%XS_hjW_ratio(j), div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(5,:) = (/ t%lhc8%XS_tthj_ratio(j), div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)

   channel_SM(1,:) = (/ t%lhc8%XS_H_SM(j) , t%BR_HZZ_SM(j) /)
   channel_SM(2,:) = (/ t%lhc8%XS_vbf_SM(j), t%BR_HZZ_SM(j) /)
   channel_SM(3,:) = (/ t%lhc8%XS_HZ_SM(j), t%BR_HZZ_SM(j) /)
   channel_SM(4,:) = (/ t%lhc8%XS_HW_SM(j), t%BR_HZZ_SM(j) /)
   channel_SM(5,:) = (/ t%lhc8%XS_ttH_SM(j), t%BR_HZZ_SM(j) /)

  case(2011026,11005,11016,11026,3478,3357,2011148,2012016) 
!   ns = 2; nb = 2; call initialise_XS_rat_BR_rat    
   nc = 4; call initialise_channel_rat_SM
        
!   XS_rat(1) = t%lhc7%XS_hj_ratio(j)   
!   XS_rat(2) = t%lhc7%XS_vbf_ratio(j)
  
!   BR_rat(1) = div( t%BR_hjWW(j) , t%BR_HWW_SM(j) ,0.0D0,1.0D0)     
!   BR_rat(2) = div(t%BR_hjZZ(j)  , t%BR_HZZ_SM(j) ,0.0D0,1.0D0)

   channel_rat(1,:) = (/ t%lhc7%XS_hj_ratio(j) , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(2,:) = (/ t%lhc7%XS_vbf_ratio(j), div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(3,:) = (/ t%lhc7%XS_hj_ratio(j) , div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(4,:) = (/ t%lhc7%XS_vbf_ratio(j), div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)

   channel_SM(1,:) = (/ t%lhc7%XS_H_SM(j) , t%BR_HWW_SM(j) /)
   channel_SM(2,:) = (/ t%lhc7%XS_vbf_SM(j), t%BR_HWW_SM(j) /)
   channel_SM(3,:) = (/ t%lhc7%XS_H_SM(j) , t%BR_HZZ_SM(j) /)
   channel_SM(4,:) = (/ t%lhc7%XS_vbf_SM(j), t%BR_HZZ_SM(j) /)

  case(11003,11014,2577,11024,1489) 
!   ns = 5; nb = 2; call initialise_XS_rat_BR_rat  
   nc = 10; call initialise_channel_rat_SM
        
!   XS_rat(1) = t%lhc7%XS_hj_ratio(j)   
!   XS_rat(2) = t%lhc7%XS_vbf_ratio(j)
!   XS_rat(3) = t%lhc7%XS_hjZ_ratio(j)  
!   XS_rat(4) = t%lhc7%XS_hjW_ratio(j) 
!   XS_rat(5) = t%lhc7%XS_tthj_ratio(j) 
  
!   BR_rat(1) = div( t%BR_hjWW(j) , t%BR_HWW_SM(j) ,0.0D0,1.0D0)
!   BR_rat(2) = div( t%BR_hjZZ(j) , t%BR_HZZ_SM(j) ,0.0D0,1.0D0)

   channel_rat(1,:) = (/ t%lhc7%XS_hj_ratio(j) , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(2,:) = (/ t%lhc7%XS_vbf_ratio(j), div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(3,:) = (/ t%lhc7%XS_hjZ_ratio(j), div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(4,:) = (/ t%lhc7%XS_hjW_ratio(j), div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(5,:) = (/ t%lhc7%XS_tthj_ratio(j), div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(6,:) = (/ t%lhc7%XS_hj_ratio(j) , div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(7,:) = (/ t%lhc7%XS_vbf_ratio(j), div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(8,:) = (/ t%lhc7%XS_hjZ_ratio(j), div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(9,:) = (/ t%lhc7%XS_hjW_ratio(j), div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(10,:) = (/ t%lhc7%XS_tthj_ratio(j), div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)

   channel_SM(1,:) = (/ t%lhc7%XS_H_SM(j) , t%BR_HWW_SM(j) /)
   channel_SM(2,:) = (/ t%lhc7%XS_vbf_SM(j), t%BR_HWW_SM(j) /)
   channel_SM(3,:) = (/ t%lhc7%XS_HZ_SM(j), t%BR_HWW_SM(j) /)
   channel_SM(4,:) = (/ t%lhc7%XS_HW_SM(j), t%BR_HWW_SM(j) /)
   channel_SM(5,:) = (/ t%lhc7%XS_ttH_SM(j), t%BR_HWW_SM(j) /)
   channel_SM(6,:) = (/ t%lhc7%XS_H_SM(j) , t%BR_HZZ_SM(j) /)
   channel_SM(7,:) = (/ t%lhc7%XS_vbf_SM(j), t%BR_HZZ_SM(j) /)
   channel_SM(8,:) = (/ t%lhc7%XS_HZ_SM(j), t%BR_HZZ_SM(j) /)
   channel_SM(9,:) = (/ t%lhc7%XS_HW_SM(j), t%BR_HZZ_SM(j) /)
   channel_SM(10,:) = (/ t%lhc7%XS_ttH_SM(j), t%BR_HZZ_SM(j) /)

  case(12042,13003)
    nc = 4; call initialise_channel_rat_SM
        
!   XS_rat(1) = t%lhc7%XS_hj_ratio(j)   
!   XS_rat(2) = t%lhc7%XS_vbf_ratio(j)
!   XS_rat(3) = t%lhc7%XS_hjZ_ratio(j)  
!   XS_rat(4) = t%lhc7%XS_hjW_ratio(j) 
!   XS_rat(5) = t%lhc7%XS_tthj_ratio(j) 
  
!   BR_rat(1) = div( t%BR_hjWW(j) , t%BR_HWW_SM(j) ,0.0D0,1.0D0)
!   BR_rat(2) = div( t%BR_hjZZ(j) , t%BR_HZZ_SM(j) ,0.0D0,1.0D0)

   channel_rat(1,:) = (/ t%lhc7%XS_hj_ratio(j) , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(2,:) = (/ t%lhc7%XS_vbf_ratio(j), div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(3,:) = (/ t%lhc7%XS_hjZ_ratio(j), div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(4,:) = (/ t%lhc7%XS_hjW_ratio(j), div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)

   channel_SM(1,:) = (/ t%lhc7%XS_H_SM(j) , t%BR_HWW_SM(j) /)
   channel_SM(2,:) = (/ t%lhc7%XS_vbf_SM(j), t%BR_HWW_SM(j) /)
   channel_SM(3,:) = (/ t%lhc7%XS_HZ_SM(j), t%BR_HWW_SM(j) /)
   channel_SM(4,:) = (/ t%lhc7%XS_HW_SM(j), t%BR_HWW_SM(j) /)

  case(2012014)
!   ns = 5; nb = 1; call initialise_XS_rat_BR_rat  
   nc = 5; call initialise_channel_rat_SM
        
!   XS_rat(1) = t%lhc7%XS_hj_ratio(j)   
!   XS_rat(2) = t%lhc7%XS_vbf_ratio(j)
!   XS_rat(3) = t%lhc7%XS_hjZ_ratio(j)  
!   XS_rat(4) = t%lhc7%XS_hjW_ratio(j) 
!   XS_rat(5) = t%lhc7%XS_tthj_ratio(j) 
  
!   BR_rat(1) = div( t%BR_hjWW(j) , t%BR_HWW_SM(j) ,0.0D0,1.0D0)
!   BR_rat(2) = div( t%BR_hjZZ(j) , t%BR_HZZ_SM(j) ,0.0D0,1.0D0)

   channel_rat(1,:) = (/ t%lhc7%XS_hj_ratio(j) , div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(2,:) = (/ t%lhc7%XS_vbf_ratio(j), div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(3,:) = (/ t%lhc7%XS_hjZ_ratio(j), div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(4,:) = (/ t%lhc7%XS_hjW_ratio(j), div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(5,:) = (/ t%lhc7%XS_tthj_ratio(j), div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)

   channel_SM(1,:) = (/ t%lhc7%XS_H_SM(j) , t%BR_Htautau_SM(j) /)
   channel_SM(2,:) = (/ t%lhc7%XS_vbf_SM(j), t%BR_Htautau_SM(j) /)
   channel_SM(3,:) = (/ t%lhc7%XS_HZ_SM(j), t%BR_Htautau_SM(j) /)
   channel_SM(4,:) = (/ t%lhc7%XS_HW_SM(j), t%BR_Htautau_SM(j) /)
   channel_SM(5,:) = (/ t%lhc7%XS_ttH_SM(j), t%BR_Htautau_SM(j) /)

  case(2012160)
!   ns = 5; nb = 1; call initialise_XS_rat_BR_rat  
   nc = 4; call initialise_channel_rat_SM
        
!   XS_rat(1) = t%lhc7%XS_hj_ratio(j)   
!   XS_rat(2) = t%lhc7%XS_vbf_ratio(j)
!   XS_rat(3) = t%lhc7%XS_hjZ_ratio(j)  
!   XS_rat(4) = t%lhc7%XS_hjW_ratio(j) 
!   XS_rat(5) = t%lhc7%XS_tthj_ratio(j) 
  
!   BR_rat(1) = div( t%BR_hjWW(j) , t%BR_HWW_SM(j) ,0.0D0,1.0D0)
!   BR_rat(2) = div( t%BR_hjZZ(j) , t%BR_HZZ_SM(j) ,0.0D0,1.0D0)

   channel_rat(1,:) = (/ t%lhc7%XS_hj_ratio(j) , div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(2,:) = (/ t%lhc7%XS_vbf_ratio(j), div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(3,:) = (/ t%lhc7%XS_hjZ_ratio(j), div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(4,:) = (/ t%lhc7%XS_hjW_ratio(j), div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)

   channel_SM(1,:) = (/ t%lhc7%XS_H_SM(j) , t%BR_Htautau_SM(j) /)
   channel_SM(2,:) = (/ t%lhc7%XS_vbf_SM(j), t%BR_Htautau_SM(j) /)
   channel_SM(3,:) = (/ t%lhc7%XS_HZ_SM(j), t%BR_Htautau_SM(j) /)
   channel_SM(4,:) = (/ t%lhc7%XS_HW_SM(j), t%BR_Htautau_SM(j) /)

  case(11009,11020,2011133) 
!   ns = 2; nb = 1; call initialise_XS_rat_BR_rat
   nc = 2; call initialise_channel_rat_SM
        
!   XS_rat(1) = t%lhc7%XS_hj_ratio(j)   
!   XS_rat(2) = t%lhc7%XS_vbf_ratio(j)
  
!   BR_rat(1) = div( t%BR_hjtautau(j) , t%BR_Htautau_SM(j) ,0.0D0,1.0D0)

   channel_rat(1,:) = (/ t%lhc7%XS_hj_ratio(j) , div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(2,:) = (/ t%lhc7%XS_vbf_ratio(j), div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)

   channel_SM(1,:) = (/ t%lhc7%XS_H_SM(j) , t%BR_Htautau_SM(j) /)
   channel_SM(2,:) = (/ t%lhc7%XS_vbf_SM(j), t%BR_Htautau_SM(j) /)

  case(110291) 
!   ns = 5; nb = 1; call initialise_XS_rat_BR_rat
   nc = 5; call initialise_channel_rat_SM
        
!   XS_rat(1) = t%lhc7%XS_hj_ratio(j)   
!   XS_rat(2) = t%lhc7%XS_vbf_ratio(j)
!   XS_rat(3) = t%lhc7%XS_hjZ_ratio(j)  
!   XS_rat(4) = t%lhc7%XS_hjW_ratio(j) 
!   XS_rat(5) = t%lhc7%XS_tthj_ratio(j) 
     
!   BR_rat(1) = div( t%BR_hjtautau(j) , t%BR_Htautau_SM(j) ,0.0D0,1.0D0)

   channel_rat(1,:) = (/ t%lhc7%XS_hj_ratio(j) , div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(2,:) = (/ t%lhc7%XS_vbf_ratio(j), div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(3,:) = (/ t%lhc7%XS_hjZ_ratio(j), div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(4,:) = (/ t%lhc7%XS_hjW_ratio(j), div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(5,:) = (/ t%lhc7%XS_tthj_ratio(j), div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)

   channel_SM(1,:) = (/ t%lhc7%XS_H_SM(j) , t%BR_Htautau_SM(j) /)
   channel_SM(2,:) = (/ t%lhc7%XS_vbf_SM(j), t%BR_Htautau_SM(j) /)
   channel_SM(3,:) = (/ t%lhc7%XS_HZ_SM(j), t%BR_Htautau_SM(j) /)
   channel_SM(4,:) = (/ t%lhc7%XS_HW_SM(j), t%BR_Htautau_SM(j) /)
   channel_SM(5,:) = (/ t%lhc7%XS_ttH_SM(j), t%BR_Htautau_SM(j) /)

  case(12043) 
   nc = 4; call initialise_channel_rat_SM
        
   channel_rat(1,:) = (/ t%lhc8%XS_hj_ratio(j) , div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(2,:) = (/ t%lhc8%XS_vbf_ratio(j), div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(3,:) = (/ t%lhc8%XS_hjZ_ratio(j), div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(4,:) = (/ t%lhc8%XS_hjW_ratio(j), div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
 
   channel_SM(1,:) = (/ t%lhc8%XS_H_SM(j) , t%BR_Htautau_SM(j) /)
   channel_SM(2,:) = (/ t%lhc8%XS_vbf_SM(j), t%BR_Htautau_SM(j) /)
   channel_SM(3,:) = (/ t%lhc8%XS_HZ_SM(j), t%BR_Htautau_SM(j) /)
   channel_SM(4,:) = (/ t%lhc8%XS_HW_SM(j), t%BR_Htautau_SM(j) /)

  case(2013010) 
   nc = 4; call initialise_channel_rat_SM
        
   channel_rat(1,:) = (/ t%lhc8%XS_hj_ratio(j) , div(t%BR_hjmumu(j),t%BR_Hmumu_SM(j),0.0D0,1.0D0) /)
   channel_rat(2,:) = (/ t%lhc8%XS_vbf_ratio(j), div(t%BR_hjmumu(j),t%BR_Hmumu_SM(j),0.0D0,1.0D0) /)
   channel_rat(3,:) = (/ t%lhc8%XS_hjZ_ratio(j), div(t%BR_hjmumu(j),t%BR_Hmumu_SM(j),0.0D0,1.0D0) /)
   channel_rat(4,:) = (/ t%lhc8%XS_hjW_ratio(j), div(t%BR_hjmumu(j),t%BR_Hmumu_SM(j),0.0D0,1.0D0) /)
 
   channel_SM(1,:) = (/ t%lhc8%XS_H_SM(j) , t%BR_Hmumu_SM(j) /)
   channel_SM(2,:) = (/ t%lhc8%XS_vbf_SM(j), t%BR_Hmumu_SM(j) /)
   channel_SM(3,:) = (/ t%lhc8%XS_HZ_SM(j), t%BR_Hmumu_SM(j) /)
   channel_SM(4,:) = (/ t%lhc8%XS_HW_SM(j), t%BR_Hmumu_SM(j) /)
 

  case(11031,2011103,11012) 
   nc = 2; call initialise_channel_rat_SM
        
   channel_rat(1,:) = (/ t%lhc7%XS_hjZ_ratio(j), div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0) /)
   channel_rat(2,:) = (/ t%lhc7%XS_hjW_ratio(j), div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0) /)

   channel_SM(1,:) = (/ t%lhc7%XS_HZ_SM(j), t%BR_Hbb_SM(j) /)
   channel_SM(2,:) = (/ t%lhc7%XS_HW_SM(j), t%BR_Hbb_SM(j) /)

  case(12044,13012) 
   nc = 2; call initialise_channel_rat_SM
        
   channel_rat(1,:) = (/ t%lhc8%XS_hjZ_ratio(j), div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0) /)
   channel_rat(2,:) = (/ t%lhc8%XS_hjW_ratio(j), div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0) /)

   channel_SM(1,:) = (/ t%lhc8%XS_HZ_SM(j), t%BR_Hbb_SM(j) /)
   channel_SM(2,:) = (/ t%lhc8%XS_HW_SM(j), t%BR_Hbb_SM(j) /)



  case(13011) 
   nc = 2; call initialise_channel_rat_SM
        
   channel_rat(1,:) = (/ t%lhc8%XS_hj_ratio(j) , div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0) /)
   channel_rat(2,:) = (/ t%lhc8%XS_vbf_ratio(j), div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0) /)

   channel_SM(1,:) = (/ t%lhc8%XS_H_SM(j) , t%BR_Hbb_SM(j) /)
   channel_SM(2,:) = (/ t%lhc8%XS_vbf_SM(j), t%BR_Hbb_SM(j) /)


  case(2012161) 
   nc = 2; call initialise_channel_rat_SM
        
   channel_rat(1,:) = (/ t%lhc8%XS_hjZ_ratio(j), div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0) /)
   channel_rat(2,:) = (/ t%lhc8%XS_hjW_ratio(j), div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0) /)

   channel_SM(1,:) = (/ t%lhc8%XS_HZ_SM(j), t%BR_Hbb_SM(j) /)
   channel_SM(2,:) = (/ t%lhc8%XS_HW_SM(j), t%BR_Hbb_SM(j) /)


  case(11021) 
!   ns = 2; nb = 1; call initialise_XS_rat_BR_rat  
   nc = 2; call initialise_channel_rat_SM
        
!   XS_rat(1) = t%lhc7%XS_hj_ratio(j)   
!   XS_rat(2) = t%lhc7%XS_vbf_ratio(j)  
 
!   BR_rat(1) = div( t%BR_hjgaga(j) , t%BR_Hgaga_SM(j) ,0.0D0,1.0D0)

   channel_rat(1,:) = (/ t%lhc7%XS_hj_ratio(j) , div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)
   channel_rat(2,:) = (/ t%lhc7%XS_vbf_ratio(j), div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)
   
   channel_SM(1,:) = (/ t%lhc7%XS_H_SM(j) , t%BR_Hgaga_SM(j) /)
   channel_SM(2,:) = (/ t%lhc7%XS_vbf_SM(j), t%BR_Hgaga_SM(j) /)
   
  case(110212)
!   ns = 3; nb = 1; call initialise_XS_rat_BR_rat
   nc = 3; call initialise_channel_rat_SM
   
!   XS_rat(1) = t%lhc7%XS_vbf_ratio(j)
!   XS_rat(2) = t%lhc7%XS_hjZ_ratio(j)  
!   XS_rat(3) = t%lhc7%XS_hjW_ratio(j) 
 
!   BR_rat(1) = div( t%BR_hjgaga(j) , t%BR_Hgaga_SM(j) ,0.0D0,1.0D0)

   channel_rat(1,:) = (/ t%lhc7%XS_vbf_ratio(j), div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)
   channel_rat(2,:) = (/ t%lhc7%XS_hjZ_ratio(j), div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)
   channel_rat(3,:) = (/ t%lhc7%XS_hjW_ratio(j), div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)

   channel_SM(1,:) = (/ t%lhc7%XS_vbf_SM(j), t%BR_Hgaga_SM(j) /)
   channel_SM(2,:) = (/ t%lhc7%XS_HZ_SM(j), t%BR_Hgaga_SM(j) /)
   channel_SM(3,:) = (/ t%lhc7%XS_HW_SM(j), t%BR_Hgaga_SM(j) /)   
   
  case(2011025,2011085,2011161,5895,1414,1487,12001,11010,11030) 
!   ns = 5; nb = 1; call initialise_XS_rat_BR_rat  
   nc = 5; call initialise_channel_rat_SM
        
!   XS_rat(1) = t%lhc7%XS_hj_ratio(j)   
!   XS_rat(2) = t%lhc7%XS_vbf_ratio(j)
!   XS_rat(3) = t%lhc7%XS_hjZ_ratio(j)  
!   XS_rat(4) = t%lhc7%XS_hjW_ratio(j) 
!   XS_rat(5) = t%lhc7%XS_tthj_ratio(j) 
  
!   BR_rat(1) = div( t%BR_hjgaga(j) , t%BR_Hgaga_SM(j) ,0.0D0,1.0D0)

   channel_rat(1,:) = (/ t%lhc7%XS_hj_ratio(j) , div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)
   channel_rat(2,:) = (/ t%lhc7%XS_vbf_ratio(j), div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)
   channel_rat(3,:) = (/ t%lhc7%XS_hjZ_ratio(j), div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)
   channel_rat(4,:) = (/ t%lhc7%XS_hjW_ratio(j), div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)
   channel_rat(5,:) = (/ t%lhc7%XS_tthj_ratio(j), div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)

   channel_SM(1,:) = (/ t%lhc7%XS_H_SM(j) , t%BR_Hgaga_SM(j) /)
   channel_SM(2,:) = (/ t%lhc7%XS_vbf_SM(j), t%BR_Hgaga_SM(j) /)
   channel_SM(3,:) = (/ t%lhc7%XS_HZ_SM(j), t%BR_Hgaga_SM(j) /)
   channel_SM(4,:) = (/ t%lhc7%XS_HW_SM(j), t%BR_Hgaga_SM(j) /)
   channel_SM(5,:) = (/ t%lhc7%XS_ttH_SM(j), t%BR_Hgaga_SM(j) /)

  case(2012091,2012168) 
!   ns = 5; nb = 1; call initialise_XS_rat_BR_rat  
   nc = 5; call initialise_channel_rat_SM
   
   channel_rat(1,:) = (/ t%lhc8%XS_hj_ratio(j) , div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)
   channel_rat(2,:) = (/ t%lhc8%XS_vbf_ratio(j), div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)
   channel_rat(3,:) = (/ t%lhc8%XS_hjZ_ratio(j), div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)
   channel_rat(4,:) = (/ t%lhc8%XS_hjW_ratio(j), div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)
   channel_rat(5,:) = (/ t%lhc8%XS_tthj_ratio(j), div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)

   channel_SM(1,:) = (/ t%lhc8%XS_H_SM(j) , t%BR_Hgaga_SM(j) /)
   channel_SM(2,:) = (/ t%lhc8%XS_vbf_SM(j), t%BR_Hgaga_SM(j) /)
   channel_SM(3,:) = (/ t%lhc8%XS_HZ_SM(j), t%BR_Hgaga_SM(j) /)
   channel_SM(4,:) = (/ t%lhc8%XS_HW_SM(j), t%BR_Hgaga_SM(j) /)
   channel_SM(5,:) = (/ t%lhc8%XS_ttH_SM(j), t%BR_Hgaga_SM(j) /)
   
      
  case(12015,13001)
    nc = 5; call initialise_channel_rat_SM
        

   channel_rat(1,:) = (/ t%lhc8%XS_hj_ratio(j) , div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)
   channel_rat(2,:) = (/ t%lhc8%XS_vbf_ratio(j), div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)
   channel_rat(3,:) = (/ t%lhc8%XS_hjZ_ratio(j), div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)
   channel_rat(4,:) = (/ t%lhc8%XS_hjW_ratio(j), div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)
   channel_rat(5,:) = (/ t%lhc8%XS_tthj_ratio(j), div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)

   channel_SM(1,:) = (/ t%lhc8%XS_H_SM(j) , t%BR_Hgaga_SM(j) /)
   channel_SM(2,:) = (/ t%lhc8%XS_vbf_SM(j), t%BR_Hgaga_SM(j) /)
   channel_SM(3,:) = (/ t%lhc8%XS_HZ_SM(j), t%BR_Hgaga_SM(j) /)
   channel_SM(4,:) = (/ t%lhc8%XS_HW_SM(j), t%BR_Hgaga_SM(j) /)
   channel_SM(5,:) = (/ t%lhc8%XS_ttH_SM(j), t%BR_Hgaga_SM(j) /)   

  case(13006,13075515,2013009)
    nc = 5; call initialise_channel_rat_SM
        
!   XS_rat(1) = t%lhc7%XS_hj_ratio(j)   
!   XS_rat(2) = t%lhc7%XS_vbf_ratio(j)
!   XS_rat(3) = t%lhc7%XS_hjZ_ratio(j)  
!   XS_rat(4) = t%lhc7%XS_hjW_ratio(j) 
!   XS_rat(5) = t%lhc7%XS_tthj_ratio(j) 
  
!   BR_rat(1) = div( t%BR_hjgaga(j) , t%BR_Hgaga_SM(j) ,0.0D0,1.0D0)

   channel_rat(1,:) = (/ t%lhc8%XS_hj_ratio(j) , div(t%BR_hjZga(j),t%BR_HZga_SM(j),0.0D0,1.0D0) /)
   channel_rat(2,:) = (/ t%lhc8%XS_vbf_ratio(j), div(t%BR_hjZga(j),t%BR_HZga_SM(j),0.0D0,1.0D0) /)
   channel_rat(3,:) = (/ t%lhc8%XS_hjZ_ratio(j), div(t%BR_hjZga(j),t%BR_HZga_SM(j),0.0D0,1.0D0) /)
   channel_rat(4,:) = (/ t%lhc8%XS_hjW_ratio(j), div(t%BR_hjZga(j),t%BR_HZga_SM(j),0.0D0,1.0D0) /)
   channel_rat(5,:) = (/ t%lhc8%XS_tthj_ratio(j), div(t%BR_hjZga(j),t%BR_HZga_SM(j),0.0D0,1.0D0) /)

   channel_SM(1,:) = (/ t%lhc8%XS_H_SM(j) , t%BR_HZga_SM(j) /)
   channel_SM(2,:) = (/ t%lhc8%XS_vbf_SM(j), t%BR_HZga_SM(j) /)
   channel_SM(3,:) = (/ t%lhc8%XS_HZ_SM(j), t%BR_HZga_SM(j) /)
   channel_SM(4,:) = (/ t%lhc8%XS_HW_SM(j), t%BR_HZga_SM(j) /)
   channel_SM(5,:) = (/ t%lhc8%XS_ttH_SM(j), t%BR_HZga_SM(j) /)   

   
  case(2748, 1408)
!   ns = 5; nb = 3; call initialise_XS_rat_BR_rat   
   nc = 15; call initialise_channel_rat_SM

!   XS_rat(1) = t%lhc7%XS_hjW_ratio(j)
!   XS_rat(2) = t%lhc7%XS_hj_ratio(j)  
!   XS_rat(3) = t%lhc7%XS_hjZ_ratio(j)  
!   XS_rat(4) = t%lhc7%XS_vbf_ratio(j)
!   XS_rat(5) = t%lhc7%XS_tthj_ratio(j)

!   BR_rat(1) = div(    t%BR_hjZZ(j) , t%BR_HZZ_SM(j)     ,0.0D0,1.0D0) 
!   BR_rat(2) = div(    t%BR_hjWW(j) , t%BR_HWW_SM(j)     ,0.0D0,1.0D0)   
!   BR_rat(3) = div(  t%BR_hjgaga(j) , t%BR_Hgaga_SM(j)   ,0.0D0,1.0D0) 

   channel_rat(1,:) = (/ t%lhc7%XS_hj_ratio(j) , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(2,:) = (/ t%lhc7%XS_vbf_ratio(j), div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(3,:) = (/ t%lhc7%XS_hjZ_ratio(j), div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(4,:) = (/ t%lhc7%XS_hjW_ratio(j), div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(5,:) = (/ t%lhc7%XS_tthj_ratio(j), div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(6,:) = (/ t%lhc7%XS_hj_ratio(j) , div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(7,:) = (/ t%lhc7%XS_vbf_ratio(j), div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(8,:) = (/ t%lhc7%XS_hjZ_ratio(j), div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(9,:) = (/ t%lhc7%XS_hjW_ratio(j), div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(10,:) = (/ t%lhc7%XS_tthj_ratio(j), div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(11,:) = (/ t%lhc7%XS_hj_ratio(j) , div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)
   channel_rat(12,:) = (/ t%lhc7%XS_vbf_ratio(j), div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)
   channel_rat(13,:) = (/ t%lhc7%XS_hjZ_ratio(j), div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)
   channel_rat(14,:) = (/ t%lhc7%XS_hjW_ratio(j), div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)
   channel_rat(15,:) = (/ t%lhc7%XS_tthj_ratio(j), div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)

   channel_SM(1,:) = (/ t%lhc7%XS_H_SM(j) , t%BR_HWW_SM(j) /)
   channel_SM(2,:) = (/ t%lhc7%XS_vbf_SM(j), t%BR_HWW_SM(j) /)
   channel_SM(3,:) = (/ t%lhc7%XS_HZ_SM(j), t%BR_HWW_SM(j) /)
   channel_SM(4,:) = (/ t%lhc7%XS_HW_SM(j), t%BR_HWW_SM(j) /)
   channel_SM(5,:) = (/ t%lhc7%XS_ttH_SM(j), t%BR_HWW_SM(j) /)
   channel_SM(6,:) = (/ t%lhc7%XS_H_SM(j) , t%BR_HZZ_SM(j) /)
   channel_SM(7,:) = (/ t%lhc7%XS_vbf_SM(j), t%BR_HZZ_SM(j) /)
   channel_SM(8,:) = (/ t%lhc7%XS_HZ_SM(j), t%BR_HZZ_SM(j) /)
   channel_SM(9,:) = (/ t%lhc7%XS_HW_SM(j), t%BR_HZZ_SM(j) /)
   channel_SM(10,:) = (/ t%lhc7%XS_ttH_SM(j), t%BR_HZZ_SM(j) /)
   channel_SM(11,:) = (/ t%lhc7%XS_H_SM(j) , t%BR_Hgaga_SM(j) /)
   channel_SM(12,:) = (/ t%lhc7%XS_vbf_SM(j), t%BR_Hgaga_SM(j) /)
   channel_SM(13,:) = (/ t%lhc7%XS_HZ_SM(j), t%BR_Hgaga_SM(j) /)
   channel_SM(14,:) = (/ t%lhc7%XS_HW_SM(j), t%BR_Hgaga_SM(j) /)
   channel_SM(15,:) = (/ t%lhc7%XS_ttH_SM(j), t%BR_Hgaga_SM(j) /)   
     
  case(11011,2011163)
!   ns = 5; nb = 4; call initialise_XS_rat_BR_rat   
   nc = 20; call initialise_channel_rat_SM

!   XS_rat(1) = t%lhc7%XS_hjW_ratio(j)
!   XS_rat(2) = t%lhc7%XS_hj_ratio(j)  
!   XS_rat(3) = t%lhc7%XS_hjZ_ratio(j)  
!   XS_rat(4) = t%lhc7%XS_vbf_ratio(j)
!   XS_rat(5) = t%lhc7%XS_tthj_ratio(j)

!   BR_rat(1) = div(    t%BR_hjZZ(j) , t%BR_HZZ_SM(j)     ,0.0D0,1.0D0) 
!   BR_rat(2) = div(    t%BR_hjWW(j) , t%BR_HWW_SM(j)     ,0.0D0,1.0D0)   
!   BR_rat(3) = div(  t%BR_hjgaga(j) , t%BR_Hgaga_SM(j)   ,0.0D0,1.0D0) 
!   BR_rat(4) = div(t%BR_hjtautau(j) , t%BR_Htautau_SM(j) ,0.0D0,1.0D0)

   channel_rat(1,:) = (/ t%lhc7%XS_hj_ratio(j) , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(2,:) = (/ t%lhc7%XS_vbf_ratio(j), div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(3,:) = (/ t%lhc7%XS_hjZ_ratio(j), div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(4,:) = (/ t%lhc7%XS_hjW_ratio(j), div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(5,:) = (/ t%lhc7%XS_tthj_ratio(j), div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(6,:) = (/ t%lhc7%XS_hj_ratio(j) , div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(7,:) = (/ t%lhc7%XS_vbf_ratio(j), div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(8,:) = (/ t%lhc7%XS_hjZ_ratio(j), div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(9,:) = (/ t%lhc7%XS_hjW_ratio(j), div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(10,:) = (/ t%lhc7%XS_tthj_ratio(j), div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(11,:) = (/ t%lhc7%XS_hj_ratio(j) , div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)
   channel_rat(12,:) = (/ t%lhc7%XS_vbf_ratio(j), div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)
   channel_rat(13,:) = (/ t%lhc7%XS_hjZ_ratio(j), div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)
   channel_rat(14,:) = (/ t%lhc7%XS_hjW_ratio(j), div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)
   channel_rat(15,:) = (/ t%lhc7%XS_tthj_ratio(j), div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)
   channel_rat(16,:) = (/ t%lhc7%XS_hj_ratio(j) , div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(17,:) = (/ t%lhc7%XS_vbf_ratio(j), div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(18,:) = (/ t%lhc7%XS_hjZ_ratio(j), div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(19,:) = (/ t%lhc7%XS_hjW_ratio(j), div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(20,:) = (/ t%lhc7%XS_tthj_ratio(j), div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)

   channel_SM(1,:) = (/ t%lhc7%XS_H_SM(j) , t%BR_HWW_SM(j) /)
   channel_SM(2,:) = (/ t%lhc7%XS_vbf_SM(j), t%BR_HWW_SM(j) /)
   channel_SM(3,:) = (/ t%lhc7%XS_HZ_SM(j), t%BR_HWW_SM(j) /)
   channel_SM(4,:) = (/ t%lhc7%XS_HW_SM(j), t%BR_HWW_SM(j) /)
   channel_SM(5,:) = (/ t%lhc7%XS_ttH_SM(j), t%BR_HWW_SM(j) /)
   channel_SM(6,:) = (/ t%lhc7%XS_H_SM(j) , t%BR_HZZ_SM(j) /)
   channel_SM(7,:) = (/ t%lhc7%XS_vbf_SM(j), t%BR_HZZ_SM(j) /)
   channel_SM(8,:) = (/ t%lhc7%XS_HZ_SM(j), t%BR_HZZ_SM(j) /)
   channel_SM(9,:) = (/ t%lhc7%XS_HW_SM(j), t%BR_HZZ_SM(j) /)
   channel_SM(10,:) = (/ t%lhc7%XS_ttH_SM(j), t%BR_HZZ_SM(j) /)
   channel_SM(11,:) = (/ t%lhc7%XS_H_SM(j) , t%BR_Hgaga_SM(j) /)
   channel_SM(12,:) = (/ t%lhc7%XS_vbf_SM(j), t%BR_Hgaga_SM(j) /)
   channel_SM(13,:) = (/ t%lhc7%XS_HZ_SM(j), t%BR_Hgaga_SM(j) /)
   channel_SM(14,:) = (/ t%lhc7%XS_HW_SM(j), t%BR_Hgaga_SM(j) /)
   channel_SM(15,:) = (/ t%lhc7%XS_ttH_SM(j), t%BR_Hgaga_SM(j) /)   
   channel_SM(16,:) = (/ t%lhc7%XS_H_SM(j) , t%BR_Htautau_SM(j) /)
   channel_SM(17,:) = (/ t%lhc7%XS_vbf_SM(j), t%BR_Htautau_SM(j) /)
   channel_SM(18,:) = (/ t%lhc7%XS_HZ_SM(j), t%BR_Htautau_SM(j) /)
   channel_SM(19,:) = (/ t%lhc7%XS_HW_SM(j), t%BR_Htautau_SM(j) /)
   channel_SM(20,:) = (/ t%lhc7%XS_ttH_SM(j), t%BR_Htautau_SM(j) /)

  case(2012012)
!   ns = 2; nb = 1; call initialise_XS_rat_BR_rat   
   nc = 2; call initialise_channel_rat_SM

   channel_rat(1,:) = (/ t%lhc7%XS_hj_ratio(j) , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(2,:) = (/ t%lhc7%XS_vbf_ratio(j), div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_SM(1,:) = (/ t%lhc7%XS_H_SM(j) , t%BR_HWW_SM(j) /)
   channel_SM(2,:) = (/ t%lhc7%XS_vbf_SM(j), t%BR_HWW_SM(j) /)

  case(2012158,2013030)
!   ns = 2; nb = 1; call initialise_XS_rat_BR_rat   
   nc = 4; call initialise_channel_rat_SM

   channel_rat(1,:) = (/ t%lhc8%XS_hj_ratio(j) , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(2,:) = (/ t%lhc8%XS_vbf_ratio(j), div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(3,:) = (/ t%lhc8%XS_hjZ_ratio(j), div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(4,:) = (/ t%lhc8%XS_hjW_ratio(j), div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_SM(1,:) = (/ t%lhc8%XS_H_SM(j) , t%BR_HWW_SM(j) /)
   channel_SM(2,:) = (/ t%lhc8%XS_vbf_SM(j), t%BR_HWW_SM(j) /)
   channel_SM(3,:) = (/ t%lhc8%XS_HZ_SM(j), t%BR_HWW_SM(j) /)
   channel_SM(4,:) = (/ t%lhc8%XS_HW_SM(j), t%BR_HWW_SM(j) /)

  case(2012135,12025)
   nc = 1; call initialise_channel_rat_SM

   channel_rat(1,:) = (/ t%lhc7%XS_tthj_ratio(j), div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0) /)
   channel_SM(1,:) = (/ t%lhc7%XS_ttH_SM(j), t%BR_Hbb_SM(j) /)  

  case(2011112)
!   ns = 5; nb = 4; call initialise_XS_rat_BR_rat   
   nc = 20; call initialise_channel_rat_SM

!   XS_rat(1) = t%lhc7%XS_hjW_ratio(j)
!   XS_rat(2) = t%lhc7%XS_hj_ratio(j)  
!   XS_rat(3) = t%lhc7%XS_hjZ_ratio(j)  
!   XS_rat(4) = t%lhc7%XS_vbf_ratio(j)
!   XS_rat(5) = t%lhc7%XS_tthj_ratio(j)

!   BR_rat(1) = div(    t%BR_hjZZ(j) , t%BR_HZZ_SM(j)     ,0.0D0,1.0D0) 
!   BR_rat(2) = div(    t%BR_hjWW(j) , t%BR_HWW_SM(j)     ,0.0D0,1.0D0)   
!   BR_rat(3) = div(  t%BR_hjgaga(j) , t%BR_Hgaga_SM(j)   ,0.0D0,1.0D0) 
!   BR_rat(4) = div(    t%BR_hjbb(j) , t%BR_Hbb_SM(j)     ,0.0D0,1.0D0)

   channel_rat(1,:) = (/ t%lhc7%XS_hj_ratio(j) , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(2,:) = (/ t%lhc7%XS_vbf_ratio(j), div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(3,:) = (/ t%lhc7%XS_hjZ_ratio(j), div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(4,:) = (/ t%lhc7%XS_hjW_ratio(j), div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(5,:) = (/ t%lhc7%XS_tthj_ratio(j), div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(6,:) = (/ t%lhc7%XS_hj_ratio(j) , div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(7,:) = (/ t%lhc7%XS_vbf_ratio(j), div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(8,:) = (/ t%lhc7%XS_hjZ_ratio(j), div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(9,:) = (/ t%lhc7%XS_hjW_ratio(j), div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(10,:) = (/ t%lhc7%XS_tthj_ratio(j), div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(11,:) = (/ t%lhc7%XS_hj_ratio(j) , div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)
   channel_rat(12,:) = (/ t%lhc7%XS_vbf_ratio(j), div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)
   channel_rat(13,:) = (/ t%lhc7%XS_hjZ_ratio(j), div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)
   channel_rat(14,:) = (/ t%lhc7%XS_hjW_ratio(j), div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)
   channel_rat(15,:) = (/ t%lhc7%XS_tthj_ratio(j), div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)
   channel_rat(16,:) = (/ t%lhc7%XS_hj_ratio(j) , div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0) /)
   channel_rat(17,:) = (/ t%lhc7%XS_vbf_ratio(j), div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0) /)
   channel_rat(18,:) = (/ t%lhc7%XS_hjZ_ratio(j), div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0) /)
   channel_rat(19,:) = (/ t%lhc7%XS_hjW_ratio(j), div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0) /)
   channel_rat(20,:) = (/ t%lhc7%XS_tthj_ratio(j), div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0) /)
   
   channel_SM(1,:) = (/ t%lhc7%XS_H_SM(j) , t%BR_HWW_SM(j) /)
   channel_SM(2,:) = (/ t%lhc7%XS_vbf_SM(j), t%BR_HWW_SM(j) /)
   channel_SM(3,:) = (/ t%lhc7%XS_HZ_SM(j), t%BR_HWW_SM(j) /)
   channel_SM(4,:) = (/ t%lhc7%XS_HW_SM(j), t%BR_HWW_SM(j) /)
   channel_SM(5,:) = (/ t%lhc7%XS_ttH_SM(j), t%BR_HWW_SM(j) /)
   channel_SM(6,:) = (/ t%lhc7%XS_H_SM(j) , t%BR_HZZ_SM(j) /)
   channel_SM(7,:) = (/ t%lhc7%XS_vbf_SM(j), t%BR_HZZ_SM(j) /)
   channel_SM(8,:) = (/ t%lhc7%XS_HZ_SM(j), t%BR_HZZ_SM(j) /)
   channel_SM(9,:) = (/ t%lhc7%XS_HW_SM(j), t%BR_HZZ_SM(j) /)
   channel_SM(10,:) = (/ t%lhc7%XS_ttH_SM(j), t%BR_HZZ_SM(j) /)
   channel_SM(11,:) = (/ t%lhc7%XS_H_SM(j) , t%BR_Hgaga_SM(j) /)
   channel_SM(12,:) = (/ t%lhc7%XS_vbf_SM(j), t%BR_Hgaga_SM(j) /)
   channel_SM(13,:) = (/ t%lhc7%XS_HZ_SM(j), t%BR_Hgaga_SM(j) /)
   channel_SM(14,:) = (/ t%lhc7%XS_HW_SM(j), t%BR_Hgaga_SM(j) /)
   channel_SM(15,:) = (/ t%lhc7%XS_ttH_SM(j), t%BR_Hgaga_SM(j) /)   
   channel_SM(16,:) = (/ t%lhc7%XS_H_SM(j) , t%BR_Hbb_SM(j) /)
   channel_SM(17,:) = (/ t%lhc7%XS_vbf_SM(j), t%BR_Hbb_SM(j) /)
   channel_SM(18,:) = (/ t%lhc7%XS_HZ_SM(j), t%BR_Hbb_SM(j) /)
   channel_SM(19,:) = (/ t%lhc7%XS_HW_SM(j), t%BR_Hbb_SM(j) /)
   channel_SM(20,:) = (/ t%lhc7%XS_ttH_SM(j), t%BR_Hbb_SM(j) /)  

  case(11022,11032,1488,12008,2011157,2012019,2011135)
!   ns = 5; nb = 5; call initialise_XS_rat_BR_rat   
   nc = 25; call initialise_channel_rat_SM

!   XS_rat(1) = t%lhc7%XS_hj_ratio(j)  
!   XS_rat(2) = t%lhc7%XS_vbf_ratio(j)
!   XS_rat(1) = t%lhc7%XS_hjW_ratio(j)
!   XS_rat(3) = t%lhc7%XS_hjZ_ratio(j)  
!   XS_rat(5) = t%lhc7%XS_tthj_ratio(j)

!   BR_rat(1) = div(    t%BR_hjZZ(j) , t%BR_HZZ_SM(j)     ,0.0D0,1.0D0) 
!   BR_rat(2) = div(    t%BR_hjWW(j) , t%BR_HWW_SM(j)     ,0.0D0,1.0D0)   
!   BR_rat(3) = div(  t%BR_hjgaga(j) , t%BR_Hgaga_SM(j)   ,0.0D0,1.0D0) 
!   BR_rat(4) = div(t%BR_hjtautau(j) , t%BR_Htautau_SM(j) ,0.0D0,1.0D0)
!   BR_rat(5) = div(    t%BR_hjbb(j) , t%BR_Hbb_SM(j)     ,0.0D0,1.0D0)   

   channel_rat(1,:) = (/ t%lhc7%XS_hj_ratio(j) , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(2,:) = (/ t%lhc7%XS_vbf_ratio(j), div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(3,:) = (/ t%lhc7%XS_hjZ_ratio(j), div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(4,:) = (/ t%lhc7%XS_hjW_ratio(j), div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(5,:) = (/ t%lhc7%XS_tthj_ratio(j), div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(6,:) = (/ t%lhc7%XS_hj_ratio(j) , div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(7,:) = (/ t%lhc7%XS_vbf_ratio(j), div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(8,:) = (/ t%lhc7%XS_hjZ_ratio(j), div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(9,:) = (/ t%lhc7%XS_hjW_ratio(j), div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(10,:) = (/ t%lhc7%XS_tthj_ratio(j), div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(11,:) = (/ t%lhc7%XS_hj_ratio(j) , div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)
   channel_rat(12,:) = (/ t%lhc7%XS_vbf_ratio(j), div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)
   channel_rat(13,:) = (/ t%lhc7%XS_hjZ_ratio(j), div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)
   channel_rat(14,:) = (/ t%lhc7%XS_hjW_ratio(j), div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)
   channel_rat(15,:) = (/ t%lhc7%XS_tthj_ratio(j), div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)
   channel_rat(16,:) = (/ t%lhc7%XS_hj_ratio(j) , div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(17,:) = (/ t%lhc7%XS_vbf_ratio(j), div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(18,:) = (/ t%lhc7%XS_hjZ_ratio(j), div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(19,:) = (/ t%lhc7%XS_hjW_ratio(j), div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(20,:) = (/ t%lhc7%XS_tthj_ratio(j), div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(21,:) = (/ t%lhc7%XS_hj_ratio(j) , div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0) /)
   channel_rat(22,:) = (/ t%lhc7%XS_vbf_ratio(j), div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0) /)
   channel_rat(23,:) = (/ t%lhc7%XS_hjZ_ratio(j), div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0) /)
   channel_rat(24,:) = (/ t%lhc7%XS_hjW_ratio(j), div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0) /)
   channel_rat(25,:) = (/ t%lhc7%XS_tthj_ratio(j), div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0) /)
   
   channel_SM(1,:) = (/ t%lhc7%XS_H_SM(j) , t%BR_HWW_SM(j) /)
   channel_SM(2,:) = (/ t%lhc7%XS_vbf_SM(j), t%BR_HWW_SM(j) /)
   channel_SM(3,:) = (/ t%lhc7%XS_HZ_SM(j), t%BR_HWW_SM(j) /)
   channel_SM(4,:) = (/ t%lhc7%XS_HW_SM(j), t%BR_HWW_SM(j) /)
   channel_SM(5,:) = (/ t%lhc7%XS_ttH_SM(j), t%BR_HWW_SM(j) /)
   channel_SM(6,:) = (/ t%lhc7%XS_H_SM(j) , t%BR_HZZ_SM(j) /)
   channel_SM(7,:) = (/ t%lhc7%XS_vbf_SM(j), t%BR_HZZ_SM(j) /)
   channel_SM(8,:) = (/ t%lhc7%XS_HZ_SM(j), t%BR_HZZ_SM(j) /)
   channel_SM(9,:) = (/ t%lhc7%XS_HW_SM(j), t%BR_HZZ_SM(j) /)
   channel_SM(10,:) = (/ t%lhc7%XS_ttH_SM(j), t%BR_HZZ_SM(j) /)
   channel_SM(11,:) = (/ t%lhc7%XS_H_SM(j) , t%BR_Hgaga_SM(j) /)
   channel_SM(12,:) = (/ t%lhc7%XS_vbf_SM(j), t%BR_Hgaga_SM(j) /)
   channel_SM(13,:) = (/ t%lhc7%XS_HZ_SM(j), t%BR_Hgaga_SM(j) /)
   channel_SM(14,:) = (/ t%lhc7%XS_HW_SM(j), t%BR_Hgaga_SM(j) /)
   channel_SM(15,:) = (/ t%lhc7%XS_ttH_SM(j), t%BR_Hgaga_SM(j) /)   
   channel_SM(16,:) = (/ t%lhc7%XS_H_SM(j) , t%BR_Htautau_SM(j) /)
   channel_SM(17,:) = (/ t%lhc7%XS_vbf_SM(j), t%BR_Htautau_SM(j) /)
   channel_SM(18,:) = (/ t%lhc7%XS_HZ_SM(j), t%BR_Htautau_SM(j) /)
   channel_SM(19,:) = (/ t%lhc7%XS_HW_SM(j), t%BR_Htautau_SM(j) /)
   channel_SM(20,:) = (/ t%lhc7%XS_ttH_SM(j), t%BR_Htautau_SM(j) /)
   channel_SM(21,:) = (/ t%lhc7%XS_H_SM(j) , t%BR_Hbb_SM(j) /)
   channel_SM(22,:) = (/ t%lhc7%XS_vbf_SM(j), t%BR_Hbb_SM(j) /)
   channel_SM(23,:) = (/ t%lhc7%XS_HZ_SM(j), t%BR_Hbb_SM(j) /)
   channel_SM(24,:) = (/ t%lhc7%XS_HW_SM(j), t%BR_Hbb_SM(j) /)
   channel_SM(25,:) = (/ t%lhc7%XS_ttH_SM(j), t%BR_Hbb_SM(j) /)   
  case(12045)
   nc = 25; call initialise_channel_rat_SM
   channel_rat(1,:) = (/ t%lhc8%XS_hj_ratio(j) , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(2,:) = (/ t%lhc8%XS_vbf_ratio(j), div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(3,:) = (/ t%lhc8%XS_hjZ_ratio(j), div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(4,:) = (/ t%lhc8%XS_hjW_ratio(j), div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(5,:) = (/ t%lhc8%XS_tthj_ratio(j), div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(6,:) = (/ t%lhc8%XS_hj_ratio(j) , div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(7,:) = (/ t%lhc8%XS_vbf_ratio(j), div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(8,:) = (/ t%lhc8%XS_hjZ_ratio(j), div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(9,:) = (/ t%lhc8%XS_hjW_ratio(j), div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(10,:) = (/ t%lhc8%XS_tthj_ratio(j), div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(11,:) = (/ t%lhc8%XS_hj_ratio(j) , div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)
   channel_rat(12,:) = (/ t%lhc8%XS_vbf_ratio(j), div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)
   channel_rat(13,:) = (/ t%lhc8%XS_hjZ_ratio(j), div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)
   channel_rat(14,:) = (/ t%lhc8%XS_hjW_ratio(j), div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)
   channel_rat(15,:) = (/ t%lhc8%XS_tthj_ratio(j), div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)
   channel_rat(16,:) = (/ t%lhc8%XS_hj_ratio(j) , div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(17,:) = (/ t%lhc8%XS_vbf_ratio(j), div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(18,:) = (/ t%lhc8%XS_hjZ_ratio(j), div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(19,:) = (/ t%lhc8%XS_hjW_ratio(j), div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(20,:) = (/ t%lhc8%XS_tthj_ratio(j), div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(21,:) = (/ t%lhc8%XS_hj_ratio(j) , div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0) /)
   channel_rat(22,:) = (/ t%lhc8%XS_vbf_ratio(j), div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0) /)
   channel_rat(23,:) = (/ t%lhc8%XS_hjZ_ratio(j), div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0) /)
   channel_rat(24,:) = (/ t%lhc8%XS_hjW_ratio(j), div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0) /)
   channel_rat(25,:) = (/ t%lhc8%XS_tthj_ratio(j), div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0) /)
   
   channel_SM(1,:) = (/ t%lhc8%XS_H_SM(j) , t%BR_HWW_SM(j) /)
   channel_SM(2,:) = (/ t%lhc8%XS_vbf_SM(j), t%BR_HWW_SM(j) /)
   channel_SM(3,:) = (/ t%lhc8%XS_HZ_SM(j), t%BR_HWW_SM(j) /)
   channel_SM(4,:) = (/ t%lhc8%XS_HW_SM(j), t%BR_HWW_SM(j) /)
   channel_SM(5,:) = (/ t%lhc8%XS_ttH_SM(j), t%BR_HWW_SM(j) /)
   channel_SM(6,:) = (/ t%lhc8%XS_H_SM(j) , t%BR_HZZ_SM(j) /)
   channel_SM(7,:) = (/ t%lhc8%XS_vbf_SM(j), t%BR_HZZ_SM(j) /)
   channel_SM(8,:) = (/ t%lhc8%XS_HZ_SM(j), t%BR_HZZ_SM(j) /)
   channel_SM(9,:) = (/ t%lhc8%XS_HW_SM(j), t%BR_HZZ_SM(j) /)
   channel_SM(10,:) = (/ t%lhc8%XS_ttH_SM(j), t%BR_HZZ_SM(j) /)
   channel_SM(11,:) = (/ t%lhc8%XS_H_SM(j) , t%BR_Hgaga_SM(j) /)
   channel_SM(12,:) = (/ t%lhc8%XS_vbf_SM(j), t%BR_Hgaga_SM(j) /)
   channel_SM(13,:) = (/ t%lhc8%XS_HZ_SM(j), t%BR_Hgaga_SM(j) /)
   channel_SM(14,:) = (/ t%lhc8%XS_HW_SM(j), t%BR_Hgaga_SM(j) /)
   channel_SM(15,:) = (/ t%lhc8%XS_ttH_SM(j), t%BR_Hgaga_SM(j) /)   
   channel_SM(16,:) = (/ t%lhc8%XS_H_SM(j) , t%BR_Htautau_SM(j) /)
   channel_SM(17,:) = (/ t%lhc8%XS_vbf_SM(j), t%BR_Htautau_SM(j) /)
   channel_SM(18,:) = (/ t%lhc8%XS_HZ_SM(j), t%BR_Htautau_SM(j) /)
   channel_SM(19,:) = (/ t%lhc8%XS_HW_SM(j), t%BR_Htautau_SM(j) /)
   channel_SM(20,:) = (/ t%lhc8%XS_ttH_SM(j), t%BR_Htautau_SM(j) /)
   channel_SM(21,:) = (/ t%lhc8%XS_H_SM(j) , t%BR_Hbb_SM(j) /)
   channel_SM(22,:) = (/ t%lhc8%XS_vbf_SM(j), t%BR_Hbb_SM(j) /)
   channel_SM(23,:) = (/ t%lhc8%XS_HZ_SM(j), t%BR_Hbb_SM(j) /)
   channel_SM(24,:) = (/ t%lhc8%XS_HW_SM(j), t%BR_Hbb_SM(j) /)
   channel_SM(25,:) = (/ t%lhc8%XS_ttH_SM(j), t%BR_Hbb_SM(j) /)    


  case(7214)
   nc = 25; call initialise_channel_rat_SM

   channel_rat(1,:) = (/ t%lhc8%XS_hj_ratio(j) , div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(2,:) = (/ t%lhc8%XS_vbf_ratio(j), div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(3,:) = (/ t%lhc8%XS_hjZ_ratio(j), div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(4,:) = (/ t%lhc8%XS_hjW_ratio(j), div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(5,:) = (/ t%lhc8%XS_tthj_ratio(j), div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0) /)
   channel_rat(6,:) = (/ t%lhc8%XS_hj_ratio(j) , div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(7,:) = (/ t%lhc8%XS_vbf_ratio(j), div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(8,:) = (/ t%lhc8%XS_hjZ_ratio(j), div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(9,:) = (/ t%lhc8%XS_hjW_ratio(j), div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(10,:) = (/ t%lhc8%XS_tthj_ratio(j), div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0) /)
   channel_rat(11,:) = (/ t%lhc8%XS_hj_ratio(j) , div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)
   channel_rat(12,:) = (/ t%lhc8%XS_vbf_ratio(j), div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)
   channel_rat(13,:) = (/ t%lhc8%XS_hjZ_ratio(j), div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)
   channel_rat(14,:) = (/ t%lhc8%XS_hjW_ratio(j), div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)
   channel_rat(15,:) = (/ t%lhc8%XS_tthj_ratio(j), div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0) /)
   channel_rat(16,:) = (/ t%lhc7%XS_hj_ratio(j) , div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(17,:) = (/ t%lhc7%XS_vbf_ratio(j), div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(18,:) = (/ t%lhc7%XS_hjZ_ratio(j), div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(19,:) = (/ t%lhc7%XS_hjW_ratio(j), div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(20,:) = (/ t%lhc7%XS_tthj_ratio(j), div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0) /)
   channel_rat(21,:) = (/ t%lhc7%XS_hj_ratio(j) , div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0) /)
   channel_rat(22,:) = (/ t%lhc7%XS_vbf_ratio(j), div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0) /)
   channel_rat(23,:) = (/ t%lhc7%XS_hjZ_ratio(j), div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0) /)
   channel_rat(24,:) = (/ t%lhc7%XS_hjW_ratio(j), div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0) /)
   channel_rat(25,:) = (/ t%lhc7%XS_tthj_ratio(j), div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0) /)
   
   channel_SM(1,:) = (/ t%lhc8%XS_H_SM(j) , t%BR_HWW_SM(j) /)
   channel_SM(2,:) = (/ t%lhc8%XS_vbf_SM(j), t%BR_HWW_SM(j) /)
   channel_SM(3,:) = (/ t%lhc8%XS_HZ_SM(j), t%BR_HWW_SM(j) /)
   channel_SM(4,:) = (/ t%lhc8%XS_HW_SM(j), t%BR_HWW_SM(j) /)
   channel_SM(5,:) = (/ t%lhc8%XS_ttH_SM(j), t%BR_HWW_SM(j) /)
   channel_SM(6,:) = (/ t%lhc8%XS_H_SM(j) , t%BR_HZZ_SM(j) /)
   channel_SM(7,:) = (/ t%lhc8%XS_vbf_SM(j), t%BR_HZZ_SM(j) /)
   channel_SM(8,:) = (/ t%lhc8%XS_HZ_SM(j), t%BR_HZZ_SM(j) /)
   channel_SM(9,:) = (/ t%lhc8%XS_HW_SM(j), t%BR_HZZ_SM(j) /)
   channel_SM(10,:) = (/ t%lhc8%XS_ttH_SM(j), t%BR_HZZ_SM(j) /)
   channel_SM(11,:) = (/ t%lhc8%XS_H_SM(j) , t%BR_Hgaga_SM(j) /)
   channel_SM(12,:) = (/ t%lhc8%XS_vbf_SM(j), t%BR_Hgaga_SM(j) /)
   channel_SM(13,:) = (/ t%lhc8%XS_HZ_SM(j), t%BR_Hgaga_SM(j) /)
   channel_SM(14,:) = (/ t%lhc8%XS_HW_SM(j), t%BR_Hgaga_SM(j) /)
   channel_SM(15,:) = (/ t%lhc8%XS_ttH_SM(j), t%BR_Hgaga_SM(j) /)   
   channel_SM(16,:) = (/ t%lhc7%XS_H_SM(j) , t%BR_Htautau_SM(j) /)
   channel_SM(17,:) = (/ t%lhc7%XS_vbf_SM(j), t%BR_Htautau_SM(j) /)
   channel_SM(18,:) = (/ t%lhc7%XS_HZ_SM(j), t%BR_Htautau_SM(j) /)
   channel_SM(19,:) = (/ t%lhc7%XS_HW_SM(j), t%BR_Htautau_SM(j) /)
   channel_SM(20,:) = (/ t%lhc7%XS_ttH_SM(j), t%BR_Htautau_SM(j) /)
   channel_SM(21,:) = (/ t%lhc7%XS_H_SM(j) , t%BR_Hbb_SM(j) /)
   channel_SM(22,:) = (/ t%lhc7%XS_vbf_SM(j), t%BR_Hbb_SM(j) /)
   channel_SM(23,:) = (/ t%lhc7%XS_HZ_SM(j), t%BR_Hbb_SM(j) /)
   channel_SM(24,:) = (/ t%lhc7%XS_HW_SM(j), t%BR_Hbb_SM(j) /)
   channel_SM(25,:) = (/ t%lhc7%XS_ttH_SM(j), t%BR_Hbb_SM(j) /)   
   
  case(5739,10574)
   ns = 1; nb = 1; call initialise_XS_rat_BR_rat   
   XS_rat(1) = t%tev%XS_tthj_ratio(j)
   BR_rat(1) = div(t%BR_hjbb(j) , t%BR_Hbb_SM(j) ,0.0D0,1.0D0) 
   if(t%CP_value(j).eq.1)then ! analysis only applies if higgs is CP even 
   else
    correct_properties=.False.
   endif 

  case(711)
   ns = 1; nb = 1; call initialise_XS_rat_BR_rat
   XS_rat(1) = t%lep%XS_bbhj_ratio(j)
   BR_rat(1) = t%BR_hjbb(j) !note *not* normalised to SM
   if(t%CP_value(j).eq.1)then ! analysis only applies if higgs is CP even 
   else
    correct_properties=.False.
   endif   

  case(713)
   ns = 1; nb = 1; call initialise_XS_rat_BR_rat 
   XS_rat(1) = t%lep%XS_bbhj_ratio(j)
   BR_rat(1) = t%BR_hjbb(j)!note *not* normalised to SM
   if(t%CP_value(j).eq.-1)then ! analysis only applies if higgs is CP odd 
   else
    correct_properties=.False.
   endif 

  case(721)
   ns = 1; nb = 1; call initialise_XS_rat_BR_rat
   XS_rat(1) = t%lep%XS_bbhj_ratio(j)
   BR_rat(1) = t%BR_hjtautau(j)!note *not* normalised to SM
   if(t%CP_value(j).eq.1)then ! analysis only applies if higgs is CP even 
   else
    correct_properties=.False.
   endif   

  case(723)
   ns = 1; nb = 1; call initialise_XS_rat_BR_rat
   XS_rat(1) = t%lep%XS_bbhj_ratio(j)
   BR_rat(1) = t%BR_hjtautau(j)!note *not* normalised to SM
   if(t%CP_value(j).eq.-1)then ! analysis only applies if higgs is CP odd 
   else
    correct_properties=.False.
   endif 

  case(731)
   ns = 1; nb = 1; call initialise_XS_rat_BR_rat 
   XS_rat(1) = t%lep%XS_tautauhj_ratio(j)
   BR_rat(1) = t%BR_hjtautau(j)!note *not* normalised to SM
   if(t%CP_value(j).eq.1)then ! analysis only applies if higgs is CP even 
   else
    correct_properties=.False.
   endif   

  case(733)
   ns = 1; nb = 1; call initialise_XS_rat_BR_rat
   XS_rat(1) = t%lep%XS_tautauhj_ratio(j)
   BR_rat(1) = t%BR_hjtautau(j)!note *not* normalised to SM
   if(t%CP_value(j).eq.-1)then ! analysis only applies if higgs is CP odd 
   else
    correct_properties=.False.
   endif 

  case(741)
   ns = 1; nb = 1; call initialise_XS_rat_BR_rat  
   XS_rat(1) = t%lep%XS_bbhj_ratio(j)
   BR_rat(1) = t%BR_hjtautau(j)!note *not* normalised to SM
   if(t%CP_value(j).eq.1)then ! analysis only applies if higgs is CP even 
   else
    correct_properties=.False.
   endif   

  case(743)
   ns = 1; nb = 1; call initialise_XS_rat_BR_rat
   XS_rat(1) = t%lep%XS_bbhj_ratio(j)
   BR_rat(1) = t%BR_hjtautau(j)!note *not* normalised to SM
   if(t%CP_value(j).eq.-1)then ! analysis only applies if higgs is CP odd 
   else
    correct_properties=.False.
   endif

  case(1811,2011094) 
   ns = 1; nb = 1; call initialise_XS_rat_BR_rat
   XS_rat(1) = t%BR_tHpjb(j) 
   BR_rat(1) = t%BR_Hpjcs(j)

   if(    (t%BR_tHpjb(j)+t%BR_tWpb       ).le.0.98D0)then
    correct_properties=.False.
   elseif((t%BR_Hpjcs(j)+t%BR_Hpjtaunu(j)).le.0.98D0)then
    correct_properties=.False.
   endif 

  case(1812,7712,8353,11002,11008) 
   ns = 1; nb = 1; call initialise_XS_rat_BR_rat
   XS_rat(1) = t%BR_tHpjb(j)
   BR_rat(1) = t%BR_Hpjtaunu(j)

   if(    (t%BR_tHpjb(j)+t%BR_tWpb).le.0.98D0)then
    correct_properties=.False.
   elseif((t%BR_Hpjcs(j)+t%BR_Hpjtaunu(j)).le.0.98D0)then
    correct_properties=.False.
   endif 

  case(2011138,2011151,2760,2013090) 
   ns = 1; nb = 1; call initialise_XS_rat_BR_rat
   XS_rat(1) = t%BR_tHpjb(j)
   BR_rat(1) = t%BR_Hpjtaunu(j)

   if(    (t%BR_tHpjb(j)+t%BR_tWpb).le.0.98D0)then
    correct_properties=.False.
   endif 

  case(1269,1270) 
   ns = 1; nb = 1; call initialise_XS_rat_BR_rat 
   XS_rat(1) = t%BR_tHpjb(j) 
   BR_rat(1) = t%BR_Hpjcs(j)

   if(    (t%BR_tHpjb(j)+t%BR_tWpb       ).le.0.98D0)then
    correct_properties=.False.
   elseif( t%BR_Hpjcs(j) .le.0.98D0)then
    correct_properties=.False.
   endif 

  case(6224) 
   ns = 1; nb = 1; call initialise_XS_rat_BR_rat 
   XS_rat(1) = t%tev%XS_hjb_ratio(j)*t%tev%XS_Hb_c4_SM(j)
   BR_rat(1) = 1.0D0

   if(    ( t%BR_hjtautau(j)+t%BR_hjbb(j) ).le.0.98D0)then
    correct_properties=.False.
   elseif( t%BR_hjtautau(j) .le.0.06D0)then
    correct_properties=.False.
   endif 
  case(6225) 
   ns = 1; nb = 1; call initialise_XS_rat_BR_rat 
   XS_rat(1) = t%tev%XS_hjb_ratio(j)*t%tev%XS_Hb_c4_SM(j)
   BR_rat(1) = 1.0D0

   if(    ( t%BR_hjtautau(j)+t%BR_hjbb(j) ).le.0.98D0)then
    correct_properties=.False.
   elseif( t%BR_hjtautau(j) .le.0.1D0)then
    correct_properties=.False.
   endif 
  case(6226) 
   ns = 1; nb = 1; call initialise_XS_rat_BR_rat 
   XS_rat(1) = t%tev%XS_hjb_ratio(j)*t%tev%XS_Hb_c4_SM(j)
   BR_rat(1) = 1.0D0

   if(    ( t%BR_hjtautau(j)+t%BR_hjbb(j) ).le.0.98D0)then
    correct_properties=.False.
   elseif( t%BR_hjtautau(j) .le.0.14D0)then
    correct_properties=.False.
   endif 

  !case(801,802,803,811,813,821)
   !ns = 1; nb = 1; call initialise_XS_rat_BR_rat      
  ! XS_rat(1) = 
   !BR_rat(1) =

   !if((t%BR_Hpjcb(j)+t%BR_Hpjcs(j)+t%BR_Hpjtaunu(j)).gt.0.98D0)then
   !else
   ! correct_properties=.False.
   !endif 

  case default
   write(*,*)'hello id=',id
   stop 'error in subroutine model_likeness (2)'
  end select

!----------------------------------------------------------
!--New model likeness check (TS 23/03/2012)
!----------------------------------------------------------
  if(allocated(channel_rat)) then
    if(nc.ne.ubound(channel_rat,dim=1))stop'error in subroutine model_likeness (3a)'
    if(nc.ne.ubound(channel_SM,dim=1))stop'error in subroutine model_likeness (3a)'
!   Check if the channels have been filled correctly
    do ic=1,nc    
     if(abs(channel_rat(ic,1)-unset).lt.1.0D-3)stop'error in subroutine model_likeness (4a)'      
     if(abs(channel_rat(ic,2)-unset).lt.1.0D-3)stop'error in subroutine model_likeness (4a)'         
     if(abs(channel_SM(ic,1)-unset).lt.1.0D-3)stop'error in subroutine model_likeness (4a)'         
     if(abs(channel_SM(ic,2)-unset).lt.1.0D-3)stop'error in subroutine model_likeness (4a)'         
    enddo
!---Eliminate irrelevant channels (=channels with very small SM prediction).
!---Construct mean value of the ratio for the relevant channels
    nc_rel=0    
    do ic=1,nc
     if(channel_SM(ic,1).gt.vsmall.and.channel_SM(ic,2).gt.vsmall) then
      nc_rel=nc_rel+1
      channel_SM(nc_rel,:)=channel_SM(ic,:)
      channel_rat(nc_rel,:)=channel_rat(ic,:)
     endif
    enddo
        
    if(nc_rel.gt.0) then
     nc=nc_rel
     call reallocate_channel_rat_SM
    endif

!--Evaluate the total SM rate expected for the (relevant) channels
    SMrate=0.
    do ic=1,nc
     SMrate=SMrate+channel_SM(ic,1)*channel_SM(ic,2)
    enddo
!--Evaluate the predicted signal strength modifier c of the model
    c=0. 
    do ic=1,nc
!----use a weighted average of the channel rate ratios     
     if(use_weight) then
      weight = div(channel_SM(ic,1)*channel_SM(ic,2),SMrate,0.0D0,1.0D9)
     else
      weight = 1.0D0/nc
     endif
     c=c+weight*channel_rat(ic,1)*channel_rat(ic,2)
    enddo

!--Evaluate the deviation of each channel rate ratio to the signal
!--strength modifier c
    allocate(dcbyc(nc))         
    do ic=1,nc    
     dcbyc(ic)= div((channel_rat(ic,1)*channel_rat(ic,2)-c),c,0.0D0,1.0D9)
    enddo

!--Do the model likeness test
    testSMratios= 1  !passes the SM-like ratios test 
    do ic=1,nc
!----Again, evaluate the weight of the channel
     if(use_weight) then
      weight = div(channel_SM(ic,1)*channel_SM(ic,2),SMrate,0.0D0,1.0D9)
     else
      weight = 1.0D0
     endif
!----Check if the channel fulfills the model likeness criteria     
!     print *, ic, channel_rat(ic,1), channel_rat(ic,2), abs(dcbyc(ic)*weight), dcbyc(ic), weight, c
     if(abs(dcbyc(ic)*weight).gt.eps)then         
      testSMratios= -1  !fails the SM-like ratios test
     endif
    enddo 

!--Write total ratio into s and b to return later.
    s=c
    b=1.  
      
    deallocate(channel_rat)
    deallocate(channel_SM)   
    deallocate(dcbyc)    

!-If channel_rat is not allocated, use old method:  
  else

    if(ns.ne.ubound(XS_rat,dim=1))stop'error in subroutine model_likeness (3a)'
    if(nb.ne.ubound(BR_rat,dim=1))stop'error in subroutine model_likeness (3b)'

    do is=1,ns    
     if(abs(XS_rat(is)-unset).lt.1.0D-3)stop'error in subroutine model_likeness (4a)'
    enddo
    do ib=1,nb    
     if(abs(BR_rat(ib)-unset).lt.1.0D-3)stop'error in subroutine model_likeness (4b)'
    enddo
  
    s=sum(XS_rat)/ns 
    b=sum(BR_rat)/nb    

    allocate(dsbys(ns))
    do is=1,ns    
     dsbys(is)= div((XS_rat(is) -s),s, 0.0D0,1.0D9)  
    enddo
  
    allocate(dbbyb(nb))
    do ib=1,nb     
     dbbyb(ib)= div((BR_rat(ib) -b),b, 0.0D0,1.0D9)
    enddo 
                         
    testSMratios= 1  !passes the SM-like ratios test 
    do is=1,ns
      do ib=1,nb 
        if(abs( dsbys(is)+dbbyb(ib)+dsbys(is)*dbbyb(ib) ).gt.eps )then 
         testSMratios= -1  !fails the SM-like ratios test
        endif
      enddo 
    enddo  

  deallocate(dsbys)
  deallocate(dbbyb)  
  
  deallocate(XS_rat)
  deallocate(BR_rat) 
  
  endif

  if(testSMratios.lt.0)correct_properties=.False.
  
  if(correct_properties)then
   model_like= 1   !passes the model-likeness test 
   sigmaXbr=s*b
  else
   model_like= -1  !fails the model-likeness test 
   sigmaXbr=0.0D0
  endif

  contains
  !----------------------------------------
  subroutine initialise_XS_rat_BR_rat

   allocate(XS_rat(ns))
   allocate(BR_rat(nb)) 
   XS_rat=unset
   BR_rat=unset
   
   allocate(XS_SM_temp(ns))
   allocate(BR_SM_temp(nb))
   XS_SM_temp=unset
   BR_SM_temp=unset
   
  end subroutine initialise_XS_rat_BR_rat
  !----------------------------------------
  subroutine initialise_channel_rat_SM

   allocate(channel_rat(nc,2))
   allocate(channel_SM(nc,2)) 
   channel_rat=unset
   channel_SM=unset

  end subroutine initialise_channel_rat_SM
  !----------------------------------------  
  subroutine reallocate_channel_rat_SM
   
   double precision, allocatable :: reallocate_array(:,:)
   allocate(reallocate_array(nc,2))
   
   reallocate_array(1:nc,:) = channel_rat(1:nc,:)
   deallocate(channel_rat)
   allocate(channel_rat(nc,2))
   channel_rat = reallocate_array

   reallocate_array(1:nc,:) = channel_SM(1:nc,:)
   deallocate(channel_SM)
   allocate(channel_SM(nc,2))
   channel_SM = reallocate_array

   deallocate(reallocate_array)
   
  end subroutine reallocate_channel_rat_SM
  !----------------------------------------  
           
 end subroutine model_likeness

 !*********************************************************** 
 subroutine fill_blank_ft1_dat(ft1,ft1_sep,vmasslower,vmasshigher,vmass_xmin,vmass_xmax,vmass_sep,valueoutsidetable) 
  ! don't forget to deallocate f_t1%dat
  use usefulbits, only : small
   implicit none
   integer :: ilower,ihigher
   double precision, intent(in) :: ft1_sep,vmasslower,vmasshigher,vmass_xmin,vmass_xmax,vmass_sep,valueoutsidetable
   type(table1) :: ft1
    
    if(abs(vmass_xmin-vmass_xmax).lt.small)stop'problem in f_from_t3 (4)'
    ft1%sep=ft1_sep

    ! we want f_t1%xmin to be lower  than x1lower
    if((vmasslower -vmass_xmin).ge.0.0D0)then
      ilower  = int((vmasslower -vmass_xmin)/vmass_sep)+1 
    else !off lower edge of table
      ilower  = int((vmasslower -vmass_xmin)/vmass_sep)+1-1 !-1 since int rounds up for negative numbers
    endif
    ihigher = int((vmasshigher-vmass_xmin)/vmass_sep)+2 ! we want f_t1%xmax to be higher than x1higher
    ft1%xmin        =  dble(ilower  - 1)*vmass_sep + vmass_xmin
    ft1%xmax        =  dble(ihigher - 1)*vmass_sep + vmass_xmin
    ft1%nx=nint((ft1%xmax-ft1%xmin)/ft1%sep)+1

    allocate(ft1%dat(ft1%nx,1)) 
     
    ft1%dat(:,1)=valueoutsidetable    
  end subroutine fill_blank_ft1_dat 
 !*********************************************************** 
 subroutine f_from_t1(t1,vmasslower,vmasshigher,sepmultfactor,datcomp, &
                    & f_t1,valueoutsidetable)
 ! Fills the f_t1 array with the information from a t1 array 
 !
 ! Do not forget to deallocate f_t1%dat later on 
 !*********************************************************** 
  use interpolate
  use usefulbits, only : small
  implicit none
  !--------------------------------------input
  type(table1), intent(in) :: t1
  double precision, intent(in) :: vmasslower,vmasshigher,valueoutsidetable
  double precision, intent(in) :: sepmultfactor
  integer, intent(in) :: datcomp
  !-----------------------------------output
  type(table1), intent(out) :: f_t1
  !-----------------------------------internal
  integer :: i
  double precision :: interpol
  double precision :: vmass,vmass_xmin,vmass_xmax,vmass_sep
  !-------------------------------------------

  if(vmasslower.gt.vmasshigher)then
    stop'problem in f_from_t1 (1)'
  endif

  f_t1%id          =  t1%id  
  f_t1%deltax      =  t1%deltax

  vmass_xmin      =  t1%xmin
  vmass_xmax      =  t1%xmax
  vmass_sep       =  t1%sep
  f_t1%sep        =  t1%sep*sepmultfactor

  call fill_blank_ft1_dat(f_t1,f_t1%sep,vmasslower,vmasshigher,vmass_xmin,vmass_xmax,vmass_sep,valueoutsidetable) 

  do i=1,ubound(f_t1%dat,dim=1)
      vmass = dble(i-1)*f_t1%sep + f_t1%xmin
     
      if(     vmass.lt.vmass_xmin-small )then
        f_t1%dat(i,1)=valueoutsidetable
      elseif( vmass.gt.vmass_xmax+small )then
        f_t1%dat(i,1)=valueoutsidetable
      else  
        call interpolate_tabletype1(vmass,t1,datcomp,interpol)
        f_t1%dat(i,1)=interpol
      endif
  enddo

 end subroutine f_from_t1  

 !*********************************************************** 
 subroutine f_from_t2(t2,m1_at_ref_point_1,m2_at_ref_point_1,m1_at_ref_point_2,m2_at_ref_point_2, &
                    & vmassm1orm2,vmasslower,vmasshigher,sepmultfactor,datcomp, &
                    & f_t1,valueoutsidetable)
 ! Fills the f_t1 array with the information from a t2 array along a line
 ! m2 = line_grad*m1  +  line_const
 !
 ! Do not forget to deallocate f_t1%dat later on 
 !*********************************************************** 
  use interpolate
  use usefulbits, only : small
  implicit none
  !--------------------------------------input
  type(table2), intent(in) :: t2
  double precision, intent(in) :: m1_at_ref_point_1,m2_at_ref_point_1,m1_at_ref_point_2,m2_at_ref_point_2
  double precision, intent(in) :: vmasslower,vmasshigher,valueoutsidetable
  double precision, intent(in) :: sepmultfactor
  integer, intent(in) :: datcomp,vmassm1orm2
  !-----------------------------------output
  type(table1), intent(out) :: f_t1
  !-----------------------------------internal
  type(table1) :: t1
  double precision :: line_grad,line_const
  integer :: i
  logical :: const_m1,const_m2
  integer :: const_m1_i,const_m2_j
  logical :: on_m1_gridline,on_m2_gridline
  double precision :: interpol,mass1,mass2
  double precision :: m1bit,m2bit
  double precision :: vmass,vmass_xmin,vmass_xmax,vmass_sep
  integer :: ftype_selection(1)
  !-------------------------------------------

  if(vmasslower.gt.vmasshigher)then
    stop'problem in f_from_t2 (1)'
  endif

  if(abs(m1_at_ref_point_1-m1_at_ref_point_2).lt.small)then
    const_m1=.True.
    !line_grad is not needed
    !line_const is not needed
  else
    const_m1=.False.
    line_grad =(m2_at_ref_point_1-m2_at_ref_point_2)/(m1_at_ref_point_1-m1_at_ref_point_2)
    line_const=(m1_at_ref_point_1*m2_at_ref_point_2-m1_at_ref_point_2*m2_at_ref_point_1) &
        &     /(m1_at_ref_point_1-m1_at_ref_point_2)
  endif

  if(abs(m2_at_ref_point_1-m2_at_ref_point_2).lt.small)then
    const_m2=.True.
  else
    const_m2=.False.
  endif

  f_t1%id          =  t2%id  
  f_t1%deltax      =  t2%deltax

  select case(vmassm1orm2)
  case(1)
    if(const_m1)stop'problem in f_from_t2 (3a)'
    vmass_xmin      =  t2%xmin1
    vmass_xmax      =  t2%xmax1
    vmass_sep       =  t2%sep1
    f_t1%sep        =  t2%sep1*sepmultfactor
  case(2)
    if(const_m2)stop'problem in f_from_t2 (3b)'
    vmass_xmin      =  t2%xmin2
    vmass_xmax      =  t2%xmax2
    vmass_sep       =  t2%sep2
    f_t1%sep        =  t2%sep2*sepmultfactor
  case default 
    stop'problem in f_from_t2 (3)'
  end select

  call fill_blank_ft1_dat(f_t1,f_t1%sep,vmasslower,vmasshigher,vmass_xmin,vmass_xmax,vmass_sep,valueoutsidetable) 

  on_m1_gridline=.False.
  if(const_m1)then
   const_m1_i=nint(  (m1_at_ref_point_1-t2%xmin1)   /t2%sep1)+1
   m1bit= m1_at_ref_point_1    -(dble(const_m1_i-1)*t2%sep1+t2%xmin1)/t2%sep1
   if(m1bit.lt.small)on_m1_gridline=.True.
  endif

  on_m2_gridline=.False.
  if(const_m2)then
   const_m2_j=nint(  (m2_at_ref_point_1-t2%xmin2)   /t2%sep2)+1
   m2bit= m2_at_ref_point_1    -(dble(const_m2_j-1)*t2%sep2+t2%xmin2)/t2%sep2
   if(m2bit.lt.small)on_m2_gridline=.True.
  endif

  ftype_selection(1)=datcomp

  if(    on_m1_gridline )then

    call fill_t1_from_t2(t2,2,const_m1_i,ftype_selection,t1)

    call f_from_t1(t1,vmasslower,vmasshigher,sepmultfactor,datcomp, &
                    & f_t1,valueoutsidetable)
    
    deallocate(t1%dat)

  elseif(on_m2_gridline )then

    call fill_t1_from_t2(t2,1,const_m2_j,ftype_selection,t1)

    call f_from_t1(t1,vmasslower,vmasshigher,sepmultfactor,datcomp, &
                    & f_t1,valueoutsidetable)

    deallocate(t1%dat)
  else
    do i=1,ubound(f_t1%dat,dim=1)
      vmass = dble(i-1)*f_t1%sep + f_t1%xmin

      if(t2%nx2.eq.1)then
          mass1 = vmass
          mass2 = t2%xmin2
      elseif(vmassm1orm2.eq.1)then
          mass1 = vmass
          mass2 = mass1*line_grad+line_const
      else  
          mass2 = vmass 
          if(const_m1)then
            mass1 = m1_at_ref_point_1
          else
            mass1 = (mass2 - line_const)/line_grad
          endif
      endif     

      if(     vmass.lt.vmass_xmin-small )then
        f_t1%dat(i,1)=valueoutsidetable
      elseif( vmass.gt.vmass_xmax+small )then
        f_t1%dat(i,1)=valueoutsidetable
      elseif((      t2%needs_M2_gt_2M1 ).and.(2.0D0*mass1>mass2+small))then
        f_t1%dat(i,1)=valueoutsidetable
      elseif((.not.(t2%needs_M2_gt_2M1)).and.(mass1>mass2+small).and.(t2%nx2.gt.1))then 
        f_t1%dat(i,1)=valueoutsidetable   
      else  
        call interpolate_tabletype2(mass1,mass2,t2,datcomp,interpol)
        f_t1%dat(i,1)=interpol
      endif
    enddo
  endif
 end subroutine f_from_t2 
 !******************************************************************     
 subroutine f_from_slices_t2(slices_t2,m1_at_ref_point_1,m2_at_ref_point_1,m1_at_ref_point_2,m2_at_ref_point_2,z, &
                    & vmassm1orm2,vmasslower,vmasshigher,sepmultfactor,datcomp, &
                    & f_t1,valueoutsidetable)
 !******************************************************************
 ! fill the f_t1 array with the information from a t3 array at constant sf along a line
 ! m2 = line_grad*m1  +  line_const
 ! do not forget to deallocate dat
  use S95tables_type3
  use interpolate
  use usefulbits, only : small
  implicit none
  type(table2), intent(in) :: slices_t2(2)
  type(table1) :: f_t1
  double precision, intent(in) :: m1_at_ref_point_1,m2_at_ref_point_1,m1_at_ref_point_2,m2_at_ref_point_2
  double precision, intent(in) :: z,vmasslower,vmasshigher,valueoutsidetable
  double precision, intent(in) :: sepmultfactor
  double precision :: line_grad,line_const
  integer, intent(in) :: datcomp,vmassm1orm2
  integer :: i
  logical :: const_m1,const_m2
  double precision :: interpol,mass1,mass2
  double precision :: vmass,vmass_xmin,vmass_xmax,vmass_sep
  double precision :: z_below,z_above

  if(vmasslower.gt.vmasshigher)then
    stop'problem in f_from_slices_t2 (1)'
  endif

  if(abs(m1_at_ref_point_1-m1_at_ref_point_2).lt.small)then
   const_m1=.True.
  else
   const_m1=.False.
  endif

  if(abs(m2_at_ref_point_1-m2_at_ref_point_2).lt.small)then
   const_m2=.True.
  else
   const_m2=.False.
  endif

  ! check if mass is within z range of table:
  if(    .not. ( (z .ge. slices_t2(1)%z-small).and.(z .le. slices_t2(2)%z+small) )  )then !#1! written in convoluted way to get the NaNs
    f_t1%id          =  slices_t2(1)%id  
    f_t1%deltax      =  slices_t2(1)%deltax

    if((slices_t2(1)%nx2.eq.1).or.(vmassm1orm2.eq.1))then
       if(const_m1)stop'problem in f_from_slices_t2 (1a)'
       vmass_xmin      =  slices_t2(1)%xmin1
       vmass_sep       =  slices_t2(1)%sep1
       f_t1%sep        =  slices_t2(1)%sep1*sepmultfactor
    else
       if(const_m2)stop'problem in f_from_slices_t2 (1b)'
       vmass_xmin      =  slices_t2(1)%xmin2
       vmass_sep       =  slices_t2(1)%sep2
       f_t1%sep        =  slices_t2(1)%sep2*sepmultfactor
    endif 

    call fill_blank_ft1_dat(f_t1,f_t1%sep,vmasslower,vmasshigher,vmass_xmin,vmass_xmax,vmass_sep,valueoutsidetable) 

  else                !#1
                  
   z_below=slices_t2(1)%z
   z_above=slices_t2(2)%z

   if(abs(z_below-z).lt.small)then !z is the same as z_below  !#2
    call f_from_t2(slices_t2(1),m1_at_ref_point_1,m2_at_ref_point_1,m1_at_ref_point_2,m2_at_ref_point_2, &
                    & vmassm1orm2,vmasslower,vmasshigher,sepmultfactor,1, &
                    & f_t1,valueoutsidetable)   

   elseif(abs(z_above-z).lt.small)then !z is the same as z_above    !#2

    call f_from_t2(slices_t2(2),m1_at_ref_point_1,m2_at_ref_point_1,m1_at_ref_point_2,m2_at_ref_point_2, &
                    & vmassm1orm2,vmasslower,vmasshigher,sepmultfactor,1, &
                    & f_t1,valueoutsidetable)                      

   else!#2   

    if(const_m1)then
     !line_grad is not needed
     !line_const is not needed
    else
     line_grad =(m2_at_ref_point_1-m2_at_ref_point_2)/(m1_at_ref_point_1-m1_at_ref_point_2)
     line_const=(m1_at_ref_point_1*m2_at_ref_point_2-m1_at_ref_point_2*m2_at_ref_point_1) &
         &     /(m1_at_ref_point_1-m1_at_ref_point_2)
    endif

    f_t1%id          =  slices_t2(1)%id  
    f_t1%deltax      =  slices_t2(1)%deltax

    if((slices_t2(1)%nx2.eq.1).or.(vmassm1orm2.eq.1))then
       vmass_xmin      =  slices_t2(1)%xmin1
       vmass_xmax      =  slices_t2(1)%xmax1
       vmass_sep       =  slices_t2(1)%sep1
       f_t1%sep        =  slices_t2(1)%sep1*sepmultfactor
    else
       if(const_m2)stop'problem in f_from_slices_t2 (3b)'
       vmass_xmin      =  slices_t2(1)%xmin2
       vmass_xmax      =  slices_t2(1)%xmax2
       vmass_sep       =  slices_t2(1)%sep2
       f_t1%sep        =  slices_t2(1)%sep2*sepmultfactor
    endif

    call fill_blank_ft1_dat(f_t1,f_t1%sep,vmasslower,vmasshigher,vmass_xmin,vmass_xmax,vmass_sep,valueoutsidetable) 
 
    do i=1,ubound(f_t1%dat,dim=1)
     vmass = dble(i-1)*f_t1%sep + f_t1%xmin 

     if(slices_t2(1)%nx2.eq.1)then
       mass1 = vmass
       mass2 = slices_t2(1)%xmin2
     else   
       select case(vmassm1orm2)
       case(1)
          mass1 = vmass
          mass2 = mass1*line_grad+line_const
       case(2)  
          mass2 = vmass 
          if(const_m1)then
            mass1 = m1_at_ref_point_1
          else
            mass1 = (mass2 - line_const)/line_grad
          endif
       case default
          stop'problem in f_from_slices_t2 (4b)'
       end select
     endif

     if(     vmass.lt.vmass_xmin-small )then
      f_t1%dat(i,1)=valueoutsidetable
     elseif( vmass.gt.vmass_xmax+small )then
      f_t1%dat(i,1)=valueoutsidetable
     elseif((slices_t2(1)%nx2.gt.1).and.(      slices_t2(1)%needs_M2_gt_2M1 ).and.(2.0D0*mass1>mass2+small))then
      f_t1%dat(i,1)=valueoutsidetable
     elseif((slices_t2(1)%nx2.gt.1).and.(.not.(slices_t2(1)%needs_M2_gt_2M1)).and.(mass1>mass2+small))then 
      f_t1%dat(i,1)=valueoutsidetable  
     else
      call interpolate_slices_t2(mass1,mass2,z,slices_t2,datcomp,interpol)
      f_t1%dat(i,1)=interpol
     endif
    enddo

   endif !#2

  endif !#1
 end subroutine f_from_slices_t2
 !******************************************************************     
 subroutine f_from_t3(t3,m1_at_ref_point_1,m2_at_ref_point_1,m1_at_ref_point_2,m2_at_ref_point_2,z, &
                    & vmassm1orm2,vmasslower,vmasshigher,sepmultfactor,datcomp, &
                    & f_t1,valueoutsidetable)
 !******************************************************************
 ! fill the f_t1 array with the information from a t3 array at constant sf along a line
 ! m2 = line_grad*m1  +  line_const
 ! do not forget to deallocate dat
  use S95tables_type3
  use interpolate
  use usefulbits, only : small
  implicit none
  type(table3), intent(in) :: t3
  type(table2) :: slices_t2(2)
  type(table1) :: f_t1
  double precision, intent(in) :: m1_at_ref_point_1,m2_at_ref_point_1,m1_at_ref_point_2,m2_at_ref_point_2
  double precision, intent(in) :: z,vmasslower,vmasshigher,valueoutsidetable
  double precision, intent(in) :: sepmultfactor
  integer, intent(in) :: datcomp,vmassm1orm2
  integer :: a
  logical :: const_m1,const_m2
  double precision :: vmass_xmin,vmass_xmax,vmass_sep
  integer :: ilow,c_zi(2),ftype_selection(1)
  double precision :: z_below,z_above

  if(vmasslower.gt.vmasshigher)then
    stop'problem in f_from_t3 (1)'
  endif

  if(abs(m1_at_ref_point_1-m1_at_ref_point_2).lt.small)then
   const_m1=.True.
  else
   const_m1=.False.
  endif

  if(abs(m2_at_ref_point_1-m2_at_ref_point_2).lt.small)then
   const_m2=.True.
  else
   const_m2=.False.
  endif

  ! check if mass is within z range of table:
  if(    .not. ( (z .ge. t3%zmin-small).and.(z .le. t3%zmax+small) )  )then !#1! written in convoluted way to get the NaNs
    f_t1%id          =  t3%id  
    f_t1%deltax      =  t3%deltax

    if((t3%nx2.eq.1).or.(vmassm1orm2.eq.1))then
       if(const_m1)stop'problem in f_from_t3 (1a)'
       vmass_xmin      =  t3%xmin1
       vmass_sep       =  t3%sep1
       f_t1%sep        =  t3%sep1*sepmultfactor
    else
       if(const_m2)stop'problem in f_from_t3 (1b)'
       vmass_xmin      =  t3%xmin2
       vmass_sep       =  t3%sep2
       f_t1%sep        =  t3%sep2*sepmultfactor
    endif 

    call fill_blank_ft1_dat(f_t1,f_t1%sep,vmasslower,vmasshigher,vmass_xmin,vmass_xmax,vmass_sep,valueoutsidetable) 

  else                !#1
                   
   ilow=int((z-t3%zmin)/t3%zsep)+1
   z_below=dble(ilow-1)*t3%zsep+t3%zmin
   z_above=z_below+t3%zsep

   if(abs(z_below-z).lt.small)then !z is the same as z_below  !#2
    c_zi= ilow
   elseif(abs(z_above-z).lt.small)then !z is the same as z_above    !#2                 
    c_zi= ilow+1          
   else !#2
    c_zi(1)= ilow
    c_zi(2)= ilow+1
   endif !#2

   ftype_selection(1)=datcomp
   call fill_slices_t2_from_slices_of_t3(t3,c_zi,ftype_selection,slices_t2)

   call f_from_slices_t2(slices_t2,m1_at_ref_point_1,m2_at_ref_point_1,m1_at_ref_point_2,m2_at_ref_point_2,z, &
                    & vmassm1orm2,vmasslower,vmasshigher,sepmultfactor,datcomp, &
                    & f_t1,valueoutsidetable)

   do a=1,2
    deallocate(slices_t2(a)%dat)
   enddo

  endif !#1
 end subroutine f_from_t3
 !************************************************************
 subroutine convolve_chisq_with_gaussian(t1,datcomp,sigma,mass,result)
 !************************************************************
 ! intergrate exp(-t1%dat(xi,1)/2)*exp(-(massx-mass)^2/(2*sigma^2))/sqrt(2*pi*sigma^2) w.r.t. x
 ! between xlower and xhigher
 ! then do -2.0D0*log to get result
 ! negative data points are invalid. They are set to zero.
  use usefulbits, only : vsmall,vvsmall,pi  !internal
  use interpolate
  use S95tables_type1
  implicit none
  type(table1),intent(in) :: t1
  integer,intent(in) :: datcomp
  double precision,intent(in) :: sigma,mass
  double precision,intent(out) :: result
  !-----------------------------------internal  
  integer :: i,ilow,ihigh,j,divisions,n,ntot
  double precision :: runningtotal,massx,datvalue,newsep
  double precision,allocatable :: newdat(:)
  double precision :: big_number_instead_of_infinity
  double precision :: dati,datiplus1
  !-------------------------------------------  
  if((datcomp.lt.lbound(t1%dat,dim=2)).or.(datcomp.gt.ubound(t1%dat,dim=2)))then
   stop'wrong datcomp inputted to subroutine convolve_with_gaussian'
  elseif(t1%nx.le.1)then
   stop'wrong t1%nx inputted to subroutine convolve_with_gaussian (2)'
  elseif(sigma.le.vsmall)then
   stop'wrong sigma inputted to subroutine convolve_with_gaussian'
  elseif(abs(t1%sep).le.vsmall)then
   stop'wrong t1%sep inputted to subroutine convolve_with_gaussian'
  endif

  big_number_instead_of_infinity=1.0D5
  divisions=5

  !do i=1,t1%nx
  ! if(t1%dat(i,datcomp).ge.big_number_instead_of_infinity)t1%dat(i,datcomp)=1.0D20
  !enddo

  n=0
  if(minval(t1%dat(:,datcomp)).lt.1.0D4)then
     ilow  = lbound(t1%dat,dim=1)
     ihigh = ubound(t1%dat,dim=1)

     if(ilow.eq.ihigh)stop'problem in subroutine convolve_with_gaussian'

     newsep=t1%sep/dble(divisions)
   
     ntot=divisions*(ihigh-ilow)+1
     allocate(newdat(ntot))
     newdat=0.0D0
     do i=ilow,ihigh
      dati=t1%dat(i,datcomp)
      if(dati.ge.0.0D0)then
         n=n+1
         massx=dble(i-1)*t1%sep+t1%xmin

         datvalue=dati
  
         newdat(n)=exp(-datvalue/2.0D0) &
                & *exp(-(massx-mass)**2.0D0/(2.0D0*sigma**2.0D0))/sqrt(2.0D0*pi*sigma**2.0D0)

         if(i.lt.ihigh)then
          datiplus1=t1%dat(i+1,datcomp)
          if(datiplus1.ge.0.0D0)then
           do j=2,divisions-1
            n=n+1
    
            massx=dble(i-1)*t1%sep+t1%xmin + dble(j-1)*newsep
            !do a=1,ihigh
            !  write(*,*)a,dble(a-1)*t1%sep+t1%xmin ,t1%dat(a,datcomp)
            !enddo
            datvalue=dati +((datiplus1-dati)/t1%sep)*dble(j-1)*newsep

  
            if(datvalue.lt.0.0D0)then !these are invalid point or places outside range of table
               datvalue=0.0D0
            endif
  
            newdat(n)=exp(-datvalue/2.0D0) &
                   & *exp(-(massx-mass)**2.0D0/(2.0D0*sigma**2.0D0))/sqrt(2.0D0*pi*sigma**2.0D0)

           enddo
          else
           do j=2,divisions-1
            n=n+1
           enddo
          endif
         else !negative data points are invalid
          do j=2,divisions-1
           n=n+1
          enddo
         endif

         !massx=dble(i-1)*t1%sep+t1%xmin
         !newdat(i)=exp(-t1%dat(i,datcomp)/2.0D0) &
         !       & *exp(-(massx-mass)**2.0D0/(2.0D0*sigma**2.0D0))/sqrt(2.0D0*pi*sigma**2.0D0)
 
      else
       do j=1,divisions-1
        n=n+1
       enddo
      endif
     enddo

     !intergrate with trapezium rule
     runningtotal=0.5D0*(newdat(1)+newdat(ntot))

     if((ntot).gt.1)then
      do n=2,ntot-1
       runningtotal=runningtotal+newdat(n)
      enddo  
     endif

     deallocate(newdat)

     if(abs(runningtotal).le.vvsmall)then
       result=  big_number_instead_of_infinity 
     else
       result= -2.0D0*log(runningtotal*newsep)
     endif

     if(result.gt.22.4D0)then !corresponds to clsb=1.0D-6, which is the lowest clsb that ppchi2 can take as input
       result=  big_number_instead_of_infinity 
     endif
   else
       result=  big_number_instead_of_infinity 
   endif
 end subroutine convolve_chisq_with_gaussian
 !************************************************************      
 function S95_t1_or_S95_t2_idfromelementnumber(ttype,tlist)
 !************************************************************ 
  implicit none
  integer :: S95_t1_or_S95_t2_idfromelementnumber
  integer,intent(in)  ::tlist
  integer,intent(in)  ::ttype

     select case(ttype)
     case(1)
       S95_t1_or_S95_t2_idfromelementnumber=S95_t1(tlist)%id
     case(2)  
       S95_t1_or_S95_t2_idfromelementnumber=S95_t2(tlist)%id
     case default
       stop'wrong input to function S95_t1_or_S95_t2_idfromelementnumber'
     end select
 end function S95_t1_or_S95_t2_idfromelementnumber
 !************************************************************      
 function S95_t1_or_S95_t2_elementnumberfromid(ttype,id)
 !************************************************************ 
  use S95tables_type1, only :t1elementnumberfromid
  use S95tables_type2, only :t2elementnumberfromid
  implicit none
  integer,intent(in)  ::id
  integer,intent(in)  ::ttype
  integer :: S95_t1_or_S95_t2_elementnumberfromid

  select case(ttype)
  case(1)
    S95_t1_or_S95_t2_elementnumberfromid= t1elementnumberfromid(S95_t1,id)
  case(2)
    S95_t1_or_S95_t2_elementnumberfromid= t2elementnumberfromid(S95_t2,id)
  case default
    stop'problem with function S95_t1_or_S95_t2_elementnumberfromid'
  end select

 end function S95_t1_or_S95_t2_elementnumberfromid
 !************************************************************      
 subroutine deallocate_S95tables
 !************************************************************
  implicit none
  !-----------------------------------internal
  integer x
  !-------------------------------------------
  do x=lbound(S95_t1,dim=1),ubound(S95_t1,dim=1)
   deallocate(S95_t1(x)%dat)
  enddo
  
  do x=lbound(S95_t2,dim=1),ubound(S95_t2,dim=1)
   deallocate(S95_t2(x)%dat)
  enddo  
  
  deallocate(S95_t1)
  deallocate(S95_t2)
  call deallocate_Exptranges
  
 end subroutine deallocate_S95tables
 !***********************************************************
 subroutine deallocate_Exptranges
 !***********************************************************
  implicit none

  if(allocated(Exptrange_Mhmin_forSMXS)) deallocate(Exptrange_Mhmin_forSMXS) 
  if(allocated(Exptrange_Mhmax_forSMXS)) deallocate(Exptrange_Mhmax_forSMXS) 
  if(allocated(Exptrange_Mhmin_forSMdecays)) deallocate(Exptrange_Mhmin_forSMdecays) 
  if(allocated(Exptrange_Mhmax_forSMdecays)) deallocate(Exptrange_Mhmax_forSMdecays) 

 end subroutine deallocate_Exptranges
 !***********************************************************

         
end module S95tables
!************************************************************
