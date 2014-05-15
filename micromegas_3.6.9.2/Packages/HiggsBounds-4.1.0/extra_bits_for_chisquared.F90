! This file is part of HiggsBounds
!  -KW
!******************************************************************
module extra_bits_for_chisquared
!******************************************************************
implicit none
!#define chatty
!#define testnewmethod
#define extrasafe

#ifdef NAGf90Fortran
 use F90_UNIX_IO, only : flush
#endif

 contains
	
 !************************************************************	
 subroutine get_chisq(theory_uncert_1s,axis_i,axis_j,sf, &
      &  y,chisq_withouttheory,chisq_withtheory)
 
  use usefulbits, only : vsmall,small, chisqcut_at_mumax
  use S95tables_type2
  use S95tables_type3
  use S95tables
  use interpolate, only : interpolate_tabletype3,  &
       &   interpolate_tabletype3_longer,interpolate_slices_t2, &
       &   interpolate_tabletype2
  double precision, intent(in) :: theory_uncert_1s
  double precision, intent(in) :: axis_i,axis_j,sf
  integer, intent(in) ::y
  double precision, intent(out) :: chisq_withouttheory,chisq_withtheory

  double precision :: sigma
  double precision :: log10sf
#ifdef chatty  
  character(LEN=10) :: createtime
  character(LEN=8) :: createdate
#endif

  type(table2) :: slices_t2(2)
  type(table2) :: mini_slices_t2(2)
  logical :: filledslice
  integer :: a,c
  integer :: ftype_selection(1) 
  double precision :: interpol_array(1)
  double precision :: massj

  double precision :: c_x(3,2) !x-values of 3 points where the value of f is known. 2nd element: above,below
  double precision :: c_y(3,2) !y-values of 3 points where the value of f is known. 2nd element: above,below
  double precision :: c_z(2) !z-values of 2 points where the value of f is known, both with same x-value and y-values
  double precision :: c_f(3,2)  !f(n,i) is f at (c_x(n,1),c_y(n,1),c_z(i))=(c_x(n,2),c_y(n,2),c_z(i)). 2nd element: above,below

  double precision :: m1_c,m2_c,log10sf_c,chisq_withouttheory_c
  double precision :: interpol
  integer :: nz
  type corner_info
    double precision :: f_new
    integer :: stat
  end type
  type(corner_info), allocatable ::corner(:)
  integer :: m,cc,zi,datcomp,ctot,i,j

  sigma=theory_uncert_1s

  if(sf.le.0.0D0)stop'wrong input to subroutine subroutine get_chisq'
  log10sf=log10(sf)

#ifdef chatty
  write(*,*)'hello axis_i,axis_j,sf',axis_i,axis_j,sf
#endif

  !---------------------------

#ifdef chatty
  call date_and_time(createdate,createtime)
  write(*,*)'beginning subroutine get_chisq    ', &
     & trim(createtime(1:2)//':'//createtime(3:4))//':'//createtime(5:6)//':'//createtime(7:8)//':'//createtime(9:10);call flush(6) 
#endif

  select case(clsb_t3(y)%type_of_assoc_table)
  case(1)
     massj=clsb_t3(y)%xmin2
  case(2)
     massj=axis_j
  case default
  end select

! The scale factor is smaller than the smallest number in the table
! This means the chi^2 is trivially zero
  if(log10sf.lt.clsb_t3(y)%zmin) then
   chisq_withouttheory = 0D0
   chisq_withtheory = 0D0
   return
  endif

!     print *, "Chisq scale factor outside table"

! 131101: This determines if the xsection is too high
! Here we should return the -999 option or set log10sf to maximal value according to the flag.
  if(log10sf.gt.clsb_t3(y)%zmax) then
   if(chisqcut_at_mumax) then
	log10sf = clsb_t3(y)%zmax
	write(*,*) 'WARNING: Signal strength too high for chisq table - using maximal value.'
   else
    chisq_withouttheory = -999.0D0
    chisq_withtheory = -999.0D0
	write(*,*) 'WARNING: Signal strength too high for chisq table - set chisq = -999.'    
	return
   endif 
  endif


  if(sigma.lt.vsmall)then

   call interpolate_tabletype3(axis_i,massj,log10sf, &
      &  clsb_t3(y),clsb_t3_dat_3rdcol_chisq,chisq_withouttheory)

   chisq_withtheory = -1.0D0

  else
   ftype_selection(1) =clsb_t3_dat_3rdcol_chisq
   call interpolate_tabletype3_longer(axis_i,massj,log10sf,clsb_t3(y),ftype_selection, &
                        &   interpol_array,slices_t2,filledslice,c_x,c_y,c_z,c_f)
   chisq_withouttheory=interpol_array(1)

#ifdef chatty
   call date_and_time(createdate,createtime)
   write(*,*)'chisq_withouttheory=',chisq_withouttheory
   write(*,*)'found c_f', &
      & trim(createtime(1:2)//':'//createtime(3:4))//':'//createtime(5:6)//':'//createtime(7:8)//':'//createtime(9:10)
   call flush(6)
   write(*,*)'c_x',c_x
   write(*,*)'c_y',c_y
   write(*,*)'c_z',c_z
   write(*,*)'c_f(1,1),c_f(2,1),c_f(3,1)',c_f(1,1),c_f(2,1),c_f(3,1)
   write(*,*)'c_f(1,2),c_f(2,2),c_f(3,2)',c_f(1,2),c_f(2,2),c_f(3,2)
   !write(*,*)'hello',ubound(slices_t2(1)%dat,dim=1),ubound(slices_t2(1)%dat,dim=2),ubound(slices_t2(1)%dat,dim=3)
   !write(*,*)'hello slices_t2(1)%dat(1,750,1)',slices_t2(1)%dat(1,750,1)
   !write(*,*)'hello slices_t2(2)%dat(1,750,1)',slices_t2(2)%dat(1,750,1)
#endif

   if(filledslice)then

    call initializetables_type2_blank(mini_slices_t2)
    do a=1,2
     mini_slices_t2(a)%xmax1       = maxval(c_x(:,a))
     mini_slices_t2(a)%xmax2       = maxval(c_y(:,a)) 
     mini_slices_t2(a)%xmin1       = minval(c_x(:,a))  
     mini_slices_t2(a)%xmin2       = minval(c_y(:,a)) 
     mini_slices_t2(a)%sep1        = slices_t2(a)%sep1
     mini_slices_t2(a)%sep2        = slices_t2(a)%sep2 
     mini_slices_t2(a)%z           = slices_t2(a)%z !only used in slices_t2
     mini_slices_t2(a)%needs_M2_gt_2M1     = slices_t2(a)%needs_M2_gt_2M1

     mini_slices_t2(a)%nx2=nint((mini_slices_t2(a)%xmax2-mini_slices_t2(a)%xmin2)/mini_slices_t2(a)%sep2)+1
     mini_slices_t2(a)%nx1=nint((mini_slices_t2(a)%xmax1-mini_slices_t2(a)%xmin1)/mini_slices_t2(a)%sep1)+1

     allocate(mini_slices_t2(a)%dat(mini_slices_t2(a)%nx2,mini_slices_t2(a)%nx1,2))    
     
    enddo
 
    if(abs(mini_slices_t2(2)%z-mini_slices_t2(1)%z).lt.small)then
     nz = 1
    else
     nz = 2
    endif
    
    ctot= mini_slices_t2(1)%nx1  *mini_slices_t2(1)%nx2 * nz
    allocate(corner(ctot))

    c=0
    do zi=1,nz
     do j=1,mini_slices_t2(zi)%nx2
      do i=1,mini_slices_t2(zi)%nx1
          c=c+1
   
          m2_c=mini_slices_t2(zi)%xmin2+dble(j-1)*mini_slices_t2(zi)%sep2
          m1_c=mini_slices_t2(zi)%xmin1+dble(i-1)*mini_slices_t2(zi)%sep1
          log10sf_c=mini_slices_t2(zi)%z

          m=0
          mini_slices_t2(zi)%dat(j,i,1)=-1.0D0

          do cc=1,3
           if((abs(c_x(cc,zi)-m1_c     )   .lt.small).and. &
            & (abs(c_y(cc,zi)-m2_c     )   .lt.small).and.&
            & (abs(c_z(zi)   -log10sf_c)   .lt.small))then
             mini_slices_t2(zi)%dat(j,i,1)=c_f(cc,zi)
             corner(c)%stat=0
             !write(*,*)'hello, c_f(cc,zi)',c_f(cc,zi),cc,zi
           endif
          enddo

          chisq_withouttheory_c=mini_slices_t2(zi)%dat(j,i,1)
          !write(*,*)'hello m1_c, m2_c, log10sf_c',m1_c, m2_c, log10sf_c

          if(corner(c)%stat.eq.0)then
             call get_chisq_with_theory_at_corners(m1_c, m2_c, log10sf_c, &
                  &   chisq_withouttheory_c, sigma, y, slices_t2, chisq_withtheory)
             corner(c)%f_new=chisq_withtheory
          else
             corner(c)%f_new=-1.0D0
          endif

          mini_slices_t2(zi)%dat(j,i,2) = corner(c)%f_new
          !write(*,*)'hello mini_slices_t2(zi)%dat(j,i,1)',mini_slices_t2(zi)%dat(j,i,1)
          !write(*,*)'hello mini_slices_t2(zi)%dat(j,i,2)',mini_slices_t2(zi)%dat(j,i,2)

      enddo
     enddo
    enddo

#ifdef chatty
      write(*,*)'mini_slices_t2(1)%dat(:,:,1)',mini_slices_t2(1)%dat(:,:,1)
      write(*,*)'~'
      write(*,*)'mini_slices_t2(2)%dat(:,:,1)',mini_slices_t2(2)%dat(:,:,1)
      write(*,*)'axis_i,massj,log10sf',axis_i,massj,log10sf
      write(*,*)'mini_slices_t2(1)%xmin1,mini_slices_t2(1)%xmin2',mini_slices_t2(1)%xmin1,mini_slices_t2(1)%xmin2
      write(*,*)'mini_slices_t2(2)%xmin1,mini_slices_t2(2)%xmin2',mini_slices_t2(2)%xmin1,mini_slices_t2(2)%xmin2
     ! call interpolate_slices_t2(axis_i,massj,log10sf,mini_slices_t2,1,interpol)
     ! write(*,*)'hello interpol',interpol
     ! call interpolate_slices_t2(axis_i,massj,log10sf,mini_slices_t2,2,interpol)
     ! write(*,*)'hello interpol',interpol
#endif

#ifdef extrasafe
    do datcomp=1,2
#else
    do datcomp=2,2
#endif

      if(nz.eq.2)then     
       call interpolate_slices_t2(axis_i,massj,log10sf,mini_slices_t2,datcomp,interpol)
      else
       call interpolate_tabletype2(axis_i,massj,mini_slices_t2(nz),datcomp,interpol)
      endif
      !write(*,*)'hello datcomp,interpol=',datcomp,interpol
#ifdef extrasafe
      if(datcomp.eq.2)then
#endif
       chisq_withtheory=interpol
#ifdef extrasafe
      else
       if(abs(chisq_withouttheory-interpol).gt.small*(abs(chisq_withouttheory)+abs(interpol)))then
         write(*,*)'hello chisq_withouttheory,interpol,nz',chisq_withouttheory,interpol,nz
         stop'problem in subroutine get_chisq (y)'
       endif
      endif
#endif
    enddo

    deallocate(corner)

    do a=1,2  
     deallocate(mini_slices_t2(a)%dat)
    enddo

   else 
     chisq_withouttheory=-1.0D0
     chisq_withtheory=-1.0D0
   endif

   do a=1,2  
    if(allocated(slices_t2(a)%dat)) deallocate(slices_t2(a)%dat)
   enddo
  endif

  if(chisq_withtheory   .gt.20.0D0)then
      chisq_withtheory   =1.0D6
  endif
  if(chisq_withouttheory.gt.20.0D0)then
      chisq_withouttheory=1.0D6
  endif
  if(chisq_withouttheory.lt.chisq_withtheory)chisq_withtheory=chisq_withouttheory

 !************************************************************ 
 end subroutine get_chisq     
 !************************************************************ 
  subroutine get_chisq_with_theory_at_corners(massi, massj, log10sf, chisq_withouttheory, sigma, y, slices_t2, chisq_withtheory)
  use usefulbits, only : vsmall
  use S95tables_type2
  use S95tables_type3
  use S95tables
  use interpolate, only : interpolate_tabletype3,interpolate_tabletype3_longer,interpolate_slices_t2

  implicit none

  double precision,intent(in) :: massi, massj, log10sf
  double precision,intent(in) :: chisq_withouttheory
  double precision,intent(in) :: sigma
  integer,intent(in) :: y   
  type(table2), intent(in) :: slices_t2(:)
  double precision,intent(out) :: chisq_withtheory 

  double precision :: Mvar
  double precision :: M0,Mexcl
#ifdef chatty
  character(LEN=10) :: createtime
  character(LEN=8) :: createdate
#endif
  double precision :: xhigher,xlower

  double precision :: sepmultfactor
  double precision :: valueoutsidetable
  integer :: vmassm1orm2
  double precision :: m1_at_ref_point_1,m2_at_ref_point_1,m1_at_ref_point_2,m2_at_ref_point_2
  double precision :: chisq_withtheory_unnormalised,result_from_convolution
  logical :: intable
  type(table1) :: f_t1_for_calc_stats

#ifdef chatty
     call date_and_time(createdate,createtime)
     write(*,*)'entered subroutine get_chisq_with_theory_at_corners         ', &
      & trim(createtime(1:2)//':'//createtime(3:4))//':'//createtime(5:6)//':'//createtime(7:8)//':'//createtime(9:10);call flush(6)
#endif

     select case(clsb_t3(y)%type_of_assoc_table)
     case(1)
        Mvar=massi
     case(2)
#ifdef testnewmethod
        Mvar=massj
#else
        Mvar=massi
#endif
     case default
     end select

      chisq_withtheory = -3.0D0

      xlower  = Mvar-sqrt(-2.0D0*sigma**2.0D0*log(1.0D-6)) !we want a range to integrate over
      xhigher = Mvar+sqrt(-2.0D0*sigma**2.0D0*log(1.0D-6)) !this is the range where 
              !(1/sqrt(2 pi sigma^2))exp(-(x-Mvar)^2/2/sigma^2) is greater than 1.0D-10*peak of gaussian

      sepmultfactor=0.05D0
      valueoutsidetable=1.0D5

#ifdef testnewmethod
   !call f_at_const_mass_and_log10sf(clsb_t3(y),xlower,xhigher,clsb_t3_dat_3rdcol_chisq,axis_i, &
   !        & axis_j,log10sf,f_t1_for_calc_stats,1.0D5,0.05D0)     
    !old method
    !call f_at_const_mass_and_log10sf(clsb_t3(elementnum), &
    !               &    clsb_t3(elementnum)%xmin2,clsb_t3(elementnum)%xmax2,  &
    !               &    datcomp,M1,M2,log10sf,f_t1, &
    !               &    valueoutsidetable,sepmultfactor)

    !new method
    !call f_from_t3(clsb_t3(elementnum),M1,clsb_t3(elementnum)%xmin2,M1,clsb_t3(elementnum)%xmax2, &
    !                & log10sf, &
    !                & vmassm1orm2,clsb_t3(elementnum)%xmin2,clsb_t3(elementnum)%xmax2,sepmultfactor,datcomp, &
    !                & f_t1,valueoutsidetable)
      select case(clsb_t3(y)%type_of_assoc_table)
      case(1)
        vmassm1orm2=1
        m1_at_ref_point_1=xlower
        m2_at_ref_point_1=massj
        m1_at_ref_point_2=xhigher
        m2_at_ref_point_2=massj
      case(2)
        vmassm1orm2=2
        m1_at_ref_point_1=massi
        m2_at_ref_point_1=xlower
        m1_at_ref_point_2=massi
        m2_at_ref_point_2=xhigher
      case default
        stop'problem in subroutine subroutine get_chisq (a)'
      end select
#else
      vmassm1orm2=1
      m1_at_ref_point_1=0.0D0
      m2_at_ref_point_1=0.0D0
      m1_at_ref_point_2=massi
      m2_at_ref_point_2=massj
    !call f_at_const_massratio_and_log10sf(clsb_t3(y),xlower,xhigher,clsb_t3_dat_3rdcol_chisq,axis_i, &
    !        & axis_j,log10sf,f_t1_for_calc_stats,0.0D0,0.05D0) 
    !call f_at_const_massratio_and_log10sf(clsb_t3(elementnum), &
    !               &    clsb_t3(elementnum)%xmin1,clsb_t3(elementnum)%xmax1,  &
    !               &    datcomp,M1,M2,log10sf,f_t1, &
    !               &    valueoutsidetable,sepmultfactor) 

    !call f_from_t3(clsb_t3(elementnum), &
    !                & 0.0D0,0.0D0, &
    !                & M1,M2, &
    !                & log10sf, &
    !                & vmassm1orm2,clsb_t3(elementnum)%xmin1,clsb_t3(elementnum)%xmax1,sepmultfactor,datcomp, &
    !                & f_t1,valueoutsidetable)
 
#endif
      call f_from_slices_t2(slices_t2,m1_at_ref_point_1,m2_at_ref_point_1,m1_at_ref_point_2,m2_at_ref_point_2, &
                    & log10sf, &
                    & vmassm1orm2,xlower,xhigher,sepmultfactor,1, &
                    & f_t1_for_calc_stats,valueoutsidetable)

#ifdef chatty
     call date_and_time(createdate,createtime)
     write(*,*)'found function. about to convolve ', &
     & trim(createtime(1:2)//':'//createtime(3:4))//':'//createtime(5:6)//':'//createtime(7:8)//':'//createtime(9:10);call flush(6)     
#endif

     call convolve_chisq_with_gaussian(f_t1_for_calc_stats,1,sigma,Mvar,result_from_convolution)

#ifdef chatty
     call date_and_time(createdate,createtime)
     write(*,*)'hello xlower,xhigher',xlower,xhigher
     !do i=1,ubound(f_t1_for_calc_stats%dat,dim=1)
     !  write(*,*)'hello2',f_t1_for_calc_stats%xmin+dble(i-1)*f_t1_for_calc_stats%sep,f_t1_for_calc_stats%dat(i,1)
     !enddo
     write(*,*)'convolved                         ', &
      & trim(createtime(1:2)//':'//createtime(3:4))//':'//createtime(5:6)//':'//createtime(7:8)//':'//createtime(9:10);call flush(6)
#endif

     chisq_withtheory_unnormalised=result_from_convolution

#ifdef chatty
       write(*,*)'chisq_withtheory_unnormalised=',chisq_withtheory_unnormalised
#endif

     if((chisq_withtheory_unnormalised.gt.20.0D0).and.(chisq_withouttheory.gt.20.0D0))then
        chisq_withtheory=1.0D6

     else
      !-------------------------------

      ! change set jjj...
      !call interpolate_tabletype3(massi,massj,log10sf, &
      !  &  clsb_t3(y),clsb_t3_dat_3rdcol_Mspecial,Mexcl) 

      !call interpolate_tabletype3(massi,massj,log10sf, &
      !  &  clsb_t3(y),clsb_t3_dat_3rdcol_M0,M0   )  

      Mexcl=gett3dat(intable,clsb_t3(y),massi,massj,log10sf,clsb_t3_dat_3rdcol_Mspecial)
      if(.not.intable)stop'problem in get_chisq_with_theory_at_corners (x1)'

      M0   =gett3dat(intable,clsb_t3(y),massi,massj,log10sf,clsb_t3_dat_3rdcol_M0)
      if(.not.intable)stop'problem in get_chisq_with_theory_at_corners (x2)'

      !... change set jjj



#ifdef chatty
   call date_and_time(createdate,createtime)
   write(*,*)'found Mexcl,M0                    ', &
      & trim(createtime(1:2)//':'//createtime(3:4))//':'//createtime(5:6)//':'//createtime(7:8)//':'//createtime(9:10)
   call flush(6)
   write(*,*)'Mexcl,M0',Mexcl,M0
#endif

      xlower  = M0-sqrt(-2.0D0*sigma**2.0D0*log(1.0D-6)) !we want a range to integrate over
      xhigher = M0+sqrt(-2.0D0*sigma**2.0D0*log(1.0D-6)) !this is the range where 
     !         !(1/sqrt(2 pi sigma^2))exp(-(x-M0)^2/2/sigma^2) is greater than 1.0D-10*peak of gaussian    

!#ifdef testnewmethod
    !call f_at_const_mass_and_log10sf(clsb_t3(y),xlower,xhigher,clsb_t3_dat_3rdcol_chisq,axis_i, &
    !     & axis_j,log10sf,f_t1_for_calc_stats,1.0D5,0.05D0) 
!#else
    !call f_at_const_massratio_and_log10sf(clsb_t3(y),xlower,xhigher,clsb_t3_dat_3rdcol_chisq,axis_i, &
    !     & axis_j,log10sf,f_t1_for_calc_stats,0.0D0,0.05D0) 
!#endif
      call f_from_slices_t2(slices_t2,m1_at_ref_point_1,m2_at_ref_point_1,m1_at_ref_point_2,m2_at_ref_point_2, &
                    & log10sf, &
                    & vmassm1orm2,xlower,xhigher,sepmultfactor,1, &
                    & f_t1_for_calc_stats,valueoutsidetable)


#ifdef chatty
      call date_and_time(createdate,createtime)
      write(*,*)'found function for normalisation  ', &
         & trim(createtime(1:2)//':'//createtime(3:4))//':'//createtime(5:6)//':'//createtime(7:8)//':'//createtime(9:10)
      call flush(6)
#endif

      call convolve_chisq_with_gaussian(f_t1_for_calc_stats,1,sigma,M0,result_from_convolution)  

#ifdef chatty
      call date_and_time(createdate,createtime)
      write(*,*)'found normalisation               ', &
         & trim(createtime(1:2)//':'//createtime(3:4))//':'//createtime(5:6)//':'//createtime(7:8)//':'//createtime(9:10)
      write(*,*)'result_from_convolution',result_from_convolution
      call flush(6)
#endif

      chisq_withtheory=chisq_withtheory_unnormalised-result_from_convolution

      if(chisq_withtheory.lt.0.0D0)then
        chisq_withtheory=0.0D0
      elseif(chisq_withtheory.gt.20.0D0)then
        chisq_withtheory=1.0D6
      endif 

      if(chisq_withouttheory.gt.2.7D0)then
        !nothing
      elseif(Mexcl.lt.vsmall)then
        !nothing
      elseif((abs(M0-Mexcl).lt.abs(Mexcl-Mvar)))then
        chisq_withtheory=0.0D0
      else
        !nothing
      endif

     endif

     if(chisq_withouttheory.lt.chisq_withtheory)chisq_withtheory=chisq_withouttheory

#ifdef chatty
     call date_and_time(createdate,createtime)
     write(*,*)'end of subroutine get_chisq_with_theory_at_corners       ', &
          & trim(createtime(1:2)//':'//createtime(3:4))//':'//createtime(5:6)//':'//createtime(7:8)//':'//createtime(9:10)
     call flush(6)
#endif

     deallocate(f_t1_for_calc_stats%dat)

  end subroutine get_chisq_with_theory_at_corners
 !************************************************************ 
end module extra_bits_for_chisquared
!******************************************************************
