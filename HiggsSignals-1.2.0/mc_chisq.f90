!--------------------------------------------------------------------
! This file is part of HiggsSignals (TS 03/03/2013)
!--------------------------------------------------------------------
module mc_chisq

use numerics
use combinatorics
use usefulbits_hs
implicit none

integer :: mc_mode = 1 
! mc_mode = 1 -> (default) Mass scan for theoretical mass uncertainty
! mc_mode = 2 ->           mu-plot smearing for theoretical mass uncertainty
contains

!--------------------------------------------------------------------
subroutine print_mc_observables
!--------------------------------------------------------------------
 implicit none
 integer :: i, j, k
 
 do i=lbound(analyses,dim=1),ubound(analyses,dim=1)
  write(*,*)
  write(*,*) '#*************************************************************************#'   
  write(*,'(A,I3,A,I3,A)') ' #                                Analysis ',i,'                             #'
  write(*,*) '#*************************************************************************#'      
  write(*,'(A25,1I10)') 'ID =', analyses(i)%table%id
  write(*,'(A25,4X,A3)') 'Collaboration =', analyses(i)%table%collaboration
  write(*,'(A25,4X,F6.2)') 'cms energy =', analyses(i)%table%energy         
  write(*,'(A25,4X,A45)') 'Reference =', analyses(i)%table%label
  write(*,'(A25,4X,A100)') 'Description =', analyses(i)%table%desc
  do j=lbound(analyses(i)%mpred%mp_Higgses,dim=1),&
  &    ubound(analyses(i)%mpred%mp_Higgses,dim=1)
   write(*,*)'#--------------------- pred. mass centered observable --------------------#'   
   write(*,'(A25,4X,F10.2)') 'Higgs mass =', analyses(i)%mpred%mp_Higgses(j)%m
   write(*,'(A25,4X,F10.2)') 'Higgs mass uncertainty =',analyses(i)%mpred%mp_Higgses(j)%dm
   write(*,'(A25,4X,F10.6)') 'Pred. signal strength =',analyses(i)%mpred%mp_Higgses(j)%mu
   write(*,'(A25,4X,I10)') 'Combined Higgses =',&
  &    size(analyses(i)%mpred%mp_Higgses(j)%Higgses,dim=1)
   do k=lbound(analyses(i)%mpred%mp_Higgses(j)%Higgses,dim=1),&
  &     ubound(analyses(i)%mpred%mp_Higgses(j)%Higgses,dim=1)
   write(*,'(A5,I2,A18,4X,3F10.4)') 'Higgs ',k,' (m, dm, mu) =', &
  &     analyses(i)%mpred%mp_Higgses(j)%Higgses(k)%m,&
  &     analyses(i)%mpred%mp_Higgses(j)%Higgses(k)%dm,&
  &     analyses(i)%mpred%mp_Higgses(j)%Higgses(k)%mu
   enddo 
   write(*,'(A25,4X,F10.6)') 'obs. mass value =', analyses(i)%mpred%mp_Higgses(j)%m_obs   
   write(*,'(A25,4X,F10.6)') 'obs. signal strength =',&
  &     analyses(i)%mpred%mp_Higgses(j)%mu_obs
   write(*,'(A25,4X,2F10.6)') 'cyan band(low,high) =',&
  & analyses(i)%mpred%mp_Higgses(j)%dmu_low_obs,analyses(i)%mpred%mp_Higgses(j)%dmu_up_obs
   write(*,'(A25,4X,2F10.6)') 'dmu0 (low,high) =',&
  &     analyses(i)%mpred%mp_Higgses(j)%dmu_low0_obs,&   
  &     analyses(i)%mpred%mp_Higgses(j)%dmu_up0_obs
   write(*,'(A25,4X,F10.6)') 'Chisq =', analyses(i)%mpred%mp_Higgses(j)%chisq   
  enddo
 enddo 

end subroutine print_mc_observables
!--------------------------------------------------------------------
subroutine print_mc_observables_essentials
!--------------------------------------------------------------------
 implicit none
 integer :: i, j, k
 
 do i=lbound(obs,dim=1),ubound(obs,dim=1)
  if(obs(i)%obstype.eq.2) then
  write(*,*)
  write(*,*) '#*************************************************************************#'   
  write(*,'(A,I3,A,I3,A)') ' #                      Mass-centered Observable ',i,&
  &'                     #'
  write(*,*) '#*************************************************************************#'      
  write(*,'(A25,1I10)') 'ID =', obs(i)%id
  write(*,'(A25,4X,A3)') 'Collaboration =', obs(i)%table%collaboration
  write(*,'(A25,4X,F6.2)') 'cms energy =', obs(i)%table%energy         
  write(*,'(A25,4X,A45)') 'Reference =', obs(i)%table%label
  write(*,'(A25,4X,A100)') 'Description =', obs(i)%table%desc
  write(*,'(A25,4X,1F5.2,A3,1F5.2)') 'luminosity, lum. error =',&
  & obs(i)%table%lumi,', ',obs(i)%table%dlumi
  write(*,'(A25,1X,A2,1F7.2,A1,1F7.2,A4,1F5.2)') 'mass range, separation =','[ ',&
  & obs(i)%table%xmin,',',obs(i)%table%xmax,' ], ',obs(i)%table%sep
  write(*,'(A25,1F10.2)') 'mass resolution =',obs(i)%table%deltam  
   write(*,*)'#------------------------ Channel information ----------------------------#'
   write(*,*)'  ID       prod.	decay 	    efficiency'
   write(*,*)'#-------------------------------------------------------------------------#'   
   do k=1, obs(i)%table%Nc
    write(*,'(1I5,5X,2A,1F15.6)') obs(i)%table%channel_id(k), &
& obs(i)%table%channel_description(k,:),obs(i)%table%channel_eff(k)
   enddo 
  write(*,*)'#-------------------------------------------------------------------------#'   
  endif
 enddo

end subroutine print_mc_observables_essentials
!--------------------------------------------------------------------
subroutine print_mc_tables_to_file
!--------------------------------------------------------------------
 use usefulbits, only : file_id_common3
 use usefulbits_hs, only : StrCompress
 implicit none
 character(LEN=100) :: formatspec
 integer :: i,kk

 formatspec='(I3,7X,I10,1X,4F7.2,1X,A3,1X,F6.2,1X,F6.2,1X,A,5X,A)'
 open(file_id_common3,file="mctables_information.txt")
 write(file_id_common3,*) "#HiggsSignals-"//trim(adjustl(HSvers))//						&
&						  " with experimental dataset '"//trim(adjustl(Exptdir))//"'" 
 write(file_id_common3,*) "#Number Analysis-ID mh_min   mh_max   mh_sep",&
 &				"   dmh_exp  collaboration   energy	luminosity   description   reference"
 write(file_id_common3,*) "#"
 kk=0
 do i=lbound(obs,dim=1),ubound(obs,dim=1)
  if(obs(i)%obstype.eq.2) then
  kk=kk+1
  write(file_id_common3,formatspec)kk,obs(i)%id,obs(i)%table%xmin,		&
  & obs(i)%table%xmax,obs(i)%table%sep,obs(i)%table%deltam,obs(i)%table%collaboration,&
  & obs(i)%table%energy,obs(i)%table%lumi,trim(strcompress(obs(i)%table%desc)),&
  & obs(i)%table%label
  endif
 enddo
 close(file_id_common3)
end subroutine print_mc_tables_to_file
!------------------------------------------------------------------------------------    
subroutine print_mc_observables_to_file
!------------------------------------------------------------------------------------    
 use usefulbits, only : file_id_common3 
 use usefulbits_hs, only : HSres
 implicit none
 character(LEN=100) :: formatspec !,formatspec2
 integer :: i,j,kk
 double precision :: mu_pull, dmu
 kk=0
 formatspec='(I3,7X,I10,1X,3F8.2,1X,I3,1X4F8.2)'
! formatspec2='(I3,7X,I10,1X,2F8.2,1X,A7,1X,A7,1X,A7,1X,5F10.4)'

 open(file_id_common3,file="mcobservables_information.txt")
 write(file_id_common3,*) "#HiggsSignals-"//trim(adjustl(HSvers))//						&
&						  " with experimental dataset '"//trim(adjustl(Exptdir))//"'" 
 write(file_id_common3,*) "#pull = (predicted - observed)/(gaussian uncertainty)"
 write(file_id_common3,*) "#Number Analysis-ID m(pred) dm(pred) mu(pred) ncomb m(obs)",&
&						  " mu(obs) dmu mu_pull"
 write(file_id_common3,*) "#"
 
 do i=lbound(analyses,dim=1),ubound(analyses,dim=1)
  do j=lbound(analyses(i)%mpred%mp_Higgses,dim=1),&
  &    ubound(analyses(i)%mpred%mp_Higgses,dim=1)
  kk=kk+1
   if(analyses(i)%mpred%mp_Higgses(j)%mu.ge.analyses(i)%mpred%mp_Higgses(j)%mu_obs) then
    dmu = analyses(i)%mpred%mp_Higgses(j)%dmu_up_obs
   else
     dmu = analyses(i)%mpred%mp_Higgses(j)%dmu_low_obs
   endif
  mu_pull = (analyses(i)%mpred%mp_Higgses(j)%mu -        &
  &          analyses(i)%mpred%mp_Higgses(j)%mu_obs)/dmu  
  write(file_id_common3,formatspec)kk,analyses(i)%id, &
&  analyses(i)%mpred%mp_Higgses(j)%m,		&
&  analyses(i)%mpred%mp_Higgses(j)%dm,       &
&  analyses(i)%mpred%mp_Higgses(j)%mu,       &
&  size(analyses(i)%mpred%mp_Higgses(j)%Higgses,dim=1), &
! combination code here!
&  analyses(i)%mpred%mp_Higgses(j)%m_obs,               &
&  analyses(i)%mpred%mp_Higgses(j)%mu_obs,              &
&  dmu,                                                 &
&  mu_pull
  enddo
 enddo 
 close(file_id_common3)
  
end subroutine print_mc_observables_to_file
!--------------------------------------------------------------------
subroutine fill_mp_obs(ii)
!--------------------------------------------------------------------
 implicit none
		
 integer, intent(in) :: ii		
 type(mp_neutHiggs), allocatable :: mp_H(:), mp_H_tmp(:)
 integer :: i,j,finished_combination,N
	 
!-Check if the Higgs is in the relevant mass region of the table
! and then tag it to be tested.  
 N=0	! = number of relevant Higgs bosons
 do i=lbound(analyses(ii)%Higgses,dim=1),ubound(analyses(ii)%Higgses,dim=1)
  if(analyses(ii)%Higgses(i)%m.ge.analyses(ii)%table%xmin.and.&
  &  analyses(ii)%Higgses(i)%m.le.analyses(ii)%table%xmax) then
   if(analyses(ii)%Higgses(i)%mp_test.eq.1) then
    analyses(ii)%Higgses(i)%mp_test=1
    N=N+1
   endif 
  else
   analyses(ii)%Higgses(i)%mp_test=0
  endif 
 enddo

!! write(*,*) "Analysis ID = ", analyses(ii)%table%id
!! write(*,*) "Higgs masses = ", analyses(ii)%Higgses(:)%m
!! write(*,*) "Higgs mass uncertainties = ", analyses(ii)%Higgses(:)%dm
!! write(*,*) "Tested with mp method = ", analyses(ii)%Higgses(:)%mp_test
!! if(N.eq.0) allocate(analyses(ii)%mpred%mp_Higgses(0))

!-First, throw all relevant single Higgses into mp_H.
 allocate(mp_H(N))
 j=0
 do i=lbound(analyses(ii)%Higgses,dim=1),ubound(analyses(ii)%Higgses,dim=1)
  if(analyses(ii)%Higgses(i)%mp_test.eq.1) then
   j=j+1
   allocate(mp_H(j)%Higgses(1))
   mp_H(j)%Higgses(1)=analyses(ii)%Higgses(i)
   mp_H(j)%m=analyses(ii)%Higgses(i)%m
   mp_H(j)%dm=analyses(ii)%Higgses(i)%dm
   mp_H(j)%mu=analyses(ii)%Higgses(i)%mu   
   mp_H(j)%mp_test=analyses(ii)%Higgses(i)%mp_test
  endif
 enddo  

!-Do the Stockholm Clustering

 finished_combination=-1
 do while(finished_combination.ne.1)    	   

!-For debugging
!  write(*,*) " "
!  write(*,*) " NEW LOOP "
!  write(*,*) " "
!  do j=lbound(mp_H,dim=1),ubound(mp_H,dim=1)
!   write(*,*) "#----------- mp_H(",j,") ------------"
!   write(*,*) "number of combined Higgses = ", size(mp_H(j)%Higgses,dim=1)
!   write(*,*) "mp_H(j)%m = ", mp_H(j)%m
!   write(*,*) "mp_H(j)%dm = ", mp_H(j)%dm
!   write(*,*) "mp_H(j)%mu = ", mp_H(j)%mu
!   write(*,*) "mp_H(j)%mp_test = ", mp_H(j)%mp_test
!  enddo 
!-END for debugging  

  call combine_Higgses(mp_H, mp_H_tmp, analyses(ii)%table, finished_combination)
  if(finished_combination.ne.1) then
   do i=lbound(mp_H,dim=1),ubound(mp_H,dim=1)
    deallocate(mp_H(i)%Higgses)
   enddo
   deallocate(mp_H)
   allocate(mp_H(size(mp_H_tmp,dim=1)))
   do i=lbound(mp_H_tmp,dim=1),ubound(mp_H_tmp,dim=1)
    allocate(mp_H(i)%Higgses(size(mp_H_tmp(i)%Higgses,dim=1)))
    mp_H(i)%Higgses = mp_H_tmp(i)%Higgses
    mp_H(i)%m = mp_H_tmp(i)%m
    mp_H(i)%dm = mp_H_tmp(i)%dm    
    mp_H(i)%mu = mp_H_tmp(i)%mu
    mp_H(i)%mp_test = mp_H_tmp(i)%mp_test
    deallocate(mp_H_tmp(i)%Higgses)
   enddo
   deallocate(mp_H_tmp) 
  endif
 enddo 

!-Two possible choices to parametrize the theoretical uncertainty: 
! (i) select the mu-obs value yielding lowest chi^2 (within box or including a chi^2
!     penalty from mass (gaussian)):
! call scan_mutable(analyses(ii)%table, mp_H)
! (ii) convolve (smear) the observed mu values with Higgs mass pdf:
 if(mc_mode.eq.2) call smear_mutable(analyses(ii)%table, mp_H)

! Fill mp_obs
 allocate(analyses(ii)%mpred%mp_Higgses(size(mp_H))) 
 do i=lbound(mp_H,dim=1),ubound(mp_H,dim=1)
  allocate(analyses(ii)%mpred%mp_Higgses(i)%Higgses(size(mp_H(i)%Higgses,dim=1)))
  analyses(ii)%mpred%mp_Higgses(i)%Higgses = mp_H(i)%Higgses
  analyses(ii)%mpred%mp_Higgses(i) = mp_H(i)
  analyses(ii)%mpred%mp_Higgses(i)%Higgscomb = Higgscomb_key(mp_H(i)%Higgses)
 enddo 

end subroutine fill_mp_obs
!----------------------------------------------------------------------
subroutine scan_mutable(table, mp_H)
!----------------------------------------------------------------------
 use datatables, only : interpolate_mutable
 use usefulbits, only : vsmall, small
 implicit none
 
 type(mutable), intent(in) :: table
 type(mp_neutHiggs), dimension(:), intent(inout) :: mp_H 

 integer :: i,j,imin,imax,col,step
 double precision :: mh, val(3), norm, boxstart, boxend, chisq, minchisq, mu_selected(3),&
 &					 mh_selected, dm, dmu

 step = 10
 call check_pdf

 imin = lbound(table%mu,dim=1)
 imax = ubound(table%mu,dim=1)

 do j=lbound(mp_H,dim=1),ubound(mp_H,dim=1)

  boxstart = mp_H(j)%m - mp_H(j)%dm
  boxend = mp_H(j)%m + mp_H(j)%dm
  if(boxstart.eq.boxend) then
!--this is the case for no theory uncertainty. Increase box bounding by one step's separation.
   boxstart = boxstart - table%sep/step
   boxend = boxend + table%sep/step           
  endif

  minchisq =vlarge
  do i= imin, step*imax  
   mh = table%xmin+(i-1)*table%sep/step
   call interpolate_mutable(val, mh, table, step)
   select case(pdf)
    case(1,3)
!---box theory uncertainty    	  
     if(mh.ge.boxstart.and.mh.le.boxend) then
      call get_lower_or_upper_value_difference(dmu,mp_H(j)%mu,val)
      chisq = (mp_H(j)%mu - val(2))**2/dmu**2
      if(chisq.le.minchisq) then
       minchisq = chisq
       mh_selected = mh
       mu_selected = val
      endif
     else
      chisq=vlarge 
     endif  
    case(2)
!---gaussian theory uncertainty
     call get_lower_or_upper_value_difference(dmu,mp_H(j)%mu,val)
     dm = mp_H(j)%dm
	 if(dm.le.small) dm = small	! Set dm to 1.0D-6 in case of no theory uncertainty.
	 chisq = (mp_H(j)%mu - val(2))**2/dmu**2 + (mp_H(j)%m - mh)**2/dm**2 
   end select 

   if(chisq.le.minchisq) then
    minchisq = chisq
    mh_selected = mh
    mu_selected = val
   endif
  enddo

!  if(norm.ne.0.) val=val/norm
!!  write(*,*) j, val, norm
  
  mp_H(j)%dmu_low_obs = abs(mu_selected(2)-mu_selected(1))
  mp_H(j)%mu_obs = mu_selected(2)
  mp_H(j)%dmu_up_obs = abs(mu_selected(3)-mu_selected(2))
  mp_H(j)%m_obs = mh_selected

!! write(*,*) "mass-centered observable at mh = ",mh_selected,"mu = ",mu_selected

 enddo

end subroutine scan_mutable
!----------------------------------------------------------------------
subroutine mass_variation_by_theory_uncertainty
!----------------------------------------------------------------------
 use datatables, only : interpolate_mutable
 use usefulbits, only : vsmall, small, np, Hneut
 implicit none

 integer :: i, ii, iii, j, jj, jjj, k, xmin, xmax, step, best_Hindex
 logical, allocatable :: testHiggs(:,:)
 logical, allocatable :: Higgstested(:)
 character(LEN=100) :: filename
 character(LEN=1) :: pdfstring

 double precision :: mh, val(3), norm, dm, dmu, sep, minchisq2,lowestmass, highestmass
 double precision, allocatable ::boxstart_max(:), boxend_min(:)
 double precision, allocatable ::boxstart(:),boxend(:),chisq(:),minchisq(:),mh_selected(:)
 
 allocate(Higgstested(np(Hneut)),chisq(np(Hneut)),boxstart(np(Hneut)),boxend(np(Hneut)), &
 &		  minchisq(np(Hneut)),mh_selected(np(Hneut)))
 step = 10

! Unphysical starting values.
 xmin=1000000.0D0
 xmax=0.0D0
 sep=1000000.0D0

 jjj=0
 do j=lbound(analyses,dim=1),ubound(analyses,dim=1)
  do jj=lbound(analyses(j)%mpred%mp_Higgses,dim=1),&
 &      ubound(analyses(j)%mpred%mp_Higgses,dim=1)
   jjj=jjj+1
   if(analyses(j)%table%xmin.le.xmin) xmin=analyses(j)%table%xmin
   if(analyses(j)%table%xmax.ge.xmax) xmax=analyses(j)%table%xmax
   if(analyses(j)%table%sep.le.sep) sep=analyses(j)%table%sep
  enddo
 enddo  
  
 allocate(testHiggs(jjj,np(Hneut)),boxstart_max(jjj),boxend_min(jjj))

! For debugging
if(additional_output) then
 write(pdfstring,'(I1)') pdf
 filename=trim(adjustl("dmth_variation_pdf_"//pdfstring//".dat"))
 open(147,file=filename)
 write(147,*) "# Mass parameter,    chi^2 at this position for each neutral Higgs boson"
 write(147,*) "#-----------------------------------------------------------------------"
endif

jjj=0
do j=lbound(analyses,dim=1),ubound(analyses,dim=1)
 do jj=lbound(analyses(j)%mpred%mp_Higgses,dim=1),&
 &     ubound(analyses(j)%mpred%mp_Higgses,dim=1)
 jjj=jjj+1
  do k=1, np(Hneut)
   testHiggs(jjj,k)=.False.
   lowestmass=vlarge
   highestmass=0.0D0
   do ii=lbound(analyses(j)%mpred%mp_Higgses(jj)%Higgses,dim=1),&
  &      ubound(analyses(j)%mpred%mp_Higgses(jj)%Higgses, dim=1)
    if(k.eq.analyses(j)%mpred%mp_Higgses(jj)%Higgses(ii)%id) testHiggs(jjj,k)=.True.
!-----In the case of combined Higgses, where there is no box overlap of the combined Higgs with
!-----the individual Higgses, we have to extend the box of the combined Higgs.
    if(analyses(j)%mpred%mp_Higgses(jj)%Higgses(ii)%m.le.lowestmass) then
     lowestmass = analyses(j)%mpred%mp_Higgses(jj)%Higgses(ii)%m
     boxstart_max(jjj) = analyses(j)%mpred%mp_Higgses(jj)%Higgses(ii)%m + &
  &                      analyses(j)%mpred%mp_Higgses(jj)%Higgses(ii)%dm - vsmall
    endif
    if(analyses(j)%mpred%mp_Higgses(jj)%Higgses(ii)%m.ge.highestmass) then
     highestmass = analyses(j)%mpred%mp_Higgses(jj)%Higgses(ii)%m
     boxend_min(jjj) = analyses(j)%mpred%mp_Higgses(jj)%Higgses(ii)%m - &
  &                    analyses(j)%mpred%mp_Higgses(jj)%Higgses(ii)%dm + vsmall 
    endif
   enddo
  enddo
 enddo
enddo   

!!write(*,*) "boxstart_max = ", boxstart_max
!!write(*,*) "boxend_min = ", boxend_min
!!write(*,*) "testHiggs = ", testHiggs

Higgstested=.False.
minchisq=vlarge-1000.0D0
! Vary mass globally for all analyses, minimize chisq.
 do i= 0, int((xmax-xmin)/(sep)*step)
  mh = xmin+i*sep/step
  chisq = 0.0D0
  jjj=0
  do j=lbound(analyses,dim=1),ubound(analyses,dim=1)
   do jj=lbound(analyses(j)%mpred%mp_Higgses,dim=1), &
 &       ubound(analyses(j)%mpred%mp_Higgses,dim=1)
    jjj=jjj+1
    call interpolate_mutable(val, mh, analyses(j)%table, step)
!---Do the variation for every neutral Higgs boson and test the relevant mpred analyses	
	do k=1, np(Hneut)
 	 if(testHiggs(jjj,k)) then
 	  Higgstested(k)=.True.
!-----n.b.: This is again based on the (possibly) combined Higgses (this is an approximation!)
!-----To do it correctly, Higgs bosons should be combined for every possible variation
      boxstart(k) = analyses(j)%mpred%mp_Higgses(jj)%m - &
 &                  analyses(j)%mpred%mp_Higgses(jj)%dm
      boxend(k)   = analyses(j)%mpred%mp_Higgses(jj)%m + &
 &                  analyses(j)%mpred%mp_Higgses(jj)%dm
      if(boxstart(k).gt.boxstart_max(jjj)) boxstart(k)=boxstart_max(jjj)
      if(boxend(k).lt.boxend_min(jjj)) boxend(k)=boxend_min(jjj)
!-----Extend box by one step size (new, has to be tested, TS 8/4/2013)
!      boxstart(k)=boxstart(k)-sep/(2*step)
!      boxend(k)=boxend(k)+sep/(2*step)      
      select case(pdf)
       case(1,3)
!------box theory uncertainty    	  
        if(mh.ge.boxstart(k).and.mh.le.boxend(k)) then
         call get_lower_or_upper_value_difference(dmu,&
  &           analyses(j)%mpred%mp_Higgses(jj)%mu,val)
         chisq(k) = chisq(k) + (analyses(j)%mpred%mp_Higgses(jj)%mu - val(2))**2/dmu**2
        else
         chisq(k)=chisq(k) + vlarge 
        endif  
       case(2)
!------gaussian theory uncertainty
        call get_lower_or_upper_value_difference(dmu,&
  &          analyses(j)%mpred%mp_Higgses(jj)%mu,val)
        dm = analyses(j)%mpred%mp_Higgses(jj)%dm
	    if(dm.le.small) dm = small	! Set dm to 1.0D-6 in case of no theory uncertainty.
	    chisq(k) = chisq(k) + (analyses(j)%mpred%mp_Higgses(jj)%mu - val(2))**2/dmu**2 + &
&				              (analyses(j)%mpred%mp_Higgses(jj)%m - mh)**2/dm**2 
       end select 
      endif 
     enddo
    enddo
   enddo  
   do k=1, np(Hneut)
    if(Higgstested(k)) then
     if(chisq(k).le.minchisq(k)) then
      minchisq(k) = chisq(k)
      mh_selected(k) = mh
!!      write(*,*) "hello, ",k, mh, chisq(k)
     endif
    else
     chisq(k)=vlarge 
    endif 
   enddo
!!-----For debugging only:
  if(additional_output) write(147,*) mh,chisq
 enddo
 if(additional_output) close(147)
!-Now that the optimal mass parameters are determined, get the observed mu values again
! and set them in the mpred observables.

  jjj=0
  do j=lbound(analyses,dim=1),ubound(analyses,dim=1)
   do jj=lbound(analyses(j)%mpred%mp_Higgses,dim=1),&
&        ubound(analyses(j)%mpred%mp_Higgses,dim=1)
    jjj=jjj+1
    minchisq2=vlarge
   	do k=1, np(Hneut)
	 testHiggs=.False.	 
	 do ii=lbound(analyses(j)%mpred%mp_Higgses(jj)%Higgses,dim=1), &
&          ubound(analyses(j)%mpred%mp_Higgses(jj)%Higgses, dim=1)
	  if(k.eq.analyses(j)%mpred%mp_Higgses(jj)%Higgses(ii)%id) testHiggs=.True.
!-Find the Higgs boson in a combination, where the chisq of the variation is lowest.
!-if no Higgs bosons are combined, we should have best_Hindex = k
	  if(minchisq(analyses(j)%mpred%mp_Higgses(jj)%Higgses(ii)%id).lt.minchisq2) then
	   minchisq2 = minchisq(analyses(j)%mpred%mp_Higgses(jj)%Higgses(ii)%id)
	   best_Hindex = analyses(j)%mpred%mp_Higgses(jj)%Higgses(ii)%id
	  endif 
	 enddo
 	 if(testHiggs(jjj,k)) then
!-One final chisq test for combined Higgses
   	  call interpolate_mutable(val, mh_selected(best_Hindex), analyses(j)%table, step)
      analyses(j)%mpred%mp_Higgses(jj)%dmu_low_obs = abs(val(2)-val(1))
      analyses(j)%mpred%mp_Higgses(jj)%mu_obs = val(2)
      analyses(j)%mpred%mp_Higgses(jj)%dmu_up_obs = abs(val(3)-val(2))
      analyses(j)%mpred%mp_Higgses(jj)%m_obs = mh_selected(best_Hindex)

!!	  write(*,*) "hello: ", val, best_Hindex, mh_selected(best_Hindex)
!n.b.: For combined Higgses, the best values obtained for the last Higgs are taken.
     endif
    enddo
   enddo
  enddo

 deallocate(Higgstested,chisq,boxstart,boxend,minchisq,mh_selected)
 deallocate(testHiggs,boxstart_max,boxend_min)

end subroutine mass_variation_by_theory_uncertainty
!----------------------------------------------------------------------
subroutine smear_mutable(table, mp_H)
!----------------------------------------------------------------------
 implicit none
 
 type(mutable), intent(in) :: table
 type(mp_neutHiggs), dimension(:), intent(inout) :: mp_H 
 
 integer :: i,j,imin,imax,col
 double precision :: mh, val(3), norm

 call check_pdf

 imin = lbound(table%mu,dim=1)
 imax = ubound(table%mu,dim=1)
  
 do j=lbound(mp_H,dim=1),ubound(mp_H,dim=1)
  val = 0.0D0
  norm = 0.0D0

  do i= imin, imax
   mh =table%xmin+(i-1)*table%sep
   do col=1,3   
    val(col) = val(col) + table%mu(i,col)*masspdf(mh, mp_H(j)%Higgses, table%sep)
   enddo
   norm = norm + masspdf(mh, mp_H(j)%Higgses, table%sep)
  enddo

  if(norm.ne.0.) val=val/norm
  
  mp_H(j)%dmu_low_obs = abs(val(2)-val(1))
  mp_H(j)%mu_obs = val(2)
  mp_H(j)%dmu_up_obs = abs(val(3)-val(2))
  mp_H(j)%m_obs = mp_H(j)%m
 enddo

end subroutine smear_mutable
!----------------------------------------------------------------------
function masspdf(mh, Higgses, sep)
! This function returns the value of the mass pdf at a given mass mh.
! In the case of combined Higgses, it is returning the sum of the
! individual mass pdfs. If the mass uncertainty is zero (delta-function)
! it returns the pdf values obtained for a mass uncertainty equal to
! the tables separation.
!----------------------------------------------------------------------
 use usefulbits, only : small
 implicit none
! 
 double precision, intent(in) :: mh, sep
 type(neutHiggs), dimension(:), intent(in) :: Higgses
 double precision :: masspdf
 integer :: i
 
 masspdf=0.0D0
 do i=lbound(Higgses,dim=1),ubound(Higgses,dim=1)
!!  write(*,*) "Mass uncertainty = ",Higgses(i)%dm
  if(Higgses(i)%dm.le.small) then
!!   masspdf =masspdf + box(mh, Higgses(i)%m, sep, sep)
   masspdf =masspdf + box(mh, Higgses(i)%m, Higgses(i)%dm, sep)
  else 
   select case (pdf)
    case(1,3) !-box, boxgaussian (==box in mp centered chi^2 method)
     masspdf = masspdf + box(mh, Higgses(i)%m, Higgses(i)%dm, sep)  
    case(2) !-gaussian
	 masspdf = masspdf + gaussian(mh, Higgses(i)%m, Higgses(i)%dm, sep)
   end select 
  endif
 enddo
 
!! write(*,*) "masspdf = ",masspdf
 
end function masspdf
!----------------------------------------------------------------------
function box(x,x0,dx,sep)
!----------------------------------------------------------------------
 use usefulbits, only : small
 implicit none
 
 double precision, intent(in) :: x, x0, dx, sep
 double precision :: box, norm
 
 norm = dx/sep
 if(norm.le.small) norm=1		! For the case of dx = 0 (delta-function)
 if((x.lt.(x0-dx-sep/2.)).or.(x.gt.(x0+dx+sep/2.))) then
  box = 0.0D0
 else
  box = 1.0D0/norm
 endif
!! write(*,*) "Norm / box = ",norm,box

  
end function box
!----------------------------------------------------------------------
function gaussian(x,x0,dx,sep)
!----------------------------------------------------------------------
 use usefulbits, only : pi
 implicit none
 
 double precision, intent(in) :: x, x0, dx, sep
 double precision :: gaussian

 gaussian=exp(-(x-x0)**2/(2*dx**2))*sep/(sqrt(2*pi)*dx)
end function
!----------------------------------------------------------------------
subroutine combine_Higgses(mp_H, mp_H_tmp, table, finished_combination)
!----------------------------------------------------------------------
 implicit none

 type(mp_neutHiggs), dimension(:), intent(inout) :: mp_H
 type(mp_neutHiggs), allocatable, intent(out) :: mp_H_tmp(:)
 type(mutable), intent(in) :: table
 integer, intent(out) :: finished_combination

 double precision :: massdiff, massdiff_tmp, m_tmp, dm_tmp, mu_tmp
 integer :: i,j, index_i, index_j


!-Determine the nearest neighboring Higgs bosons
 massdiff=100000.0D0
 do i=lbound(mp_H,dim=1),ubound(mp_H,dim=1)
  do j=lbound(mp_H,dim=1),ubound(mp_H,dim=1)
   if(i.ne.j) then
    massdiff_tmp = abs(mp_H(i)%m - mp_H(j)%m)
    if(massdiff_tmp.le.massdiff) then
     massdiff = massdiff_tmp
     index_i = i
     index_j = j
    endif
   endif
  enddo
 enddo 
 
 if(massdiff.le.table%deltam) then
  finished_combination=-1
!-Do a gaussian average of the two Higgs bosons to determine uncertainty and mass position,
!-and add the signal strength modifiers. Untag the combined Higgs bosons.
  if(mp_H(index_i)%dm.eq.0.0D0.and.mp_H(index_j)%dm.eq.0.0D0) then
   dm_tmp = 0.0D0
  else
   dm_tmp = mp_H(index_i)%dm*mp_H(index_j)%dm/									&
   &		   sqrt(mp_H(index_i)%dm**2+mp_H(index_j)%dm**2)
  endif 
  if(dm_tmp.ne.0.0D0) then
  m_tmp  = dm_tmp**2*( mp_H(index_i)%m/mp_H(index_i)%dm**2 +					&
  					   mp_H(index_j)%m/mp_H(index_j)%dm**2	)
  else
   if(mp_H(index_i)%dm.eq.0.0D0.and.mp_H(index_j)%dm.ne.0.0D0) then
    m_tmp = mp_H(index_i)%m
   else if(mp_H(index_i)%dm.ne.0.0D0.and.mp_H(index_j)%dm.eq.0.0D0) then
    m_tmp = mp_H(index_j)%m 
   else
    m_tmp = (mp_H(index_i)%m + mp_H(index_j)%m )/2.0D0
   endif
  endif 
!!  write(*,*) "m_tmp = ", m_tmp
  mu_tmp = mp_H(index_i)%mu + mp_H(index_j)%mu
  mp_H(index_i)%mp_test=0
  mp_H(index_j)%mp_test=0  
!--Fill the temporary mp_neutHiggs object first with the combination, then with the 
!--remaining Higgses
  allocate(mp_H_tmp(size(mp_H,dim=1)-1))
  mp_H_tmp(1)%m = m_tmp
  mp_H_tmp(1)%dm = dm_tmp
  mp_H_tmp(1)%mu = mu_tmp
  mp_H_tmp(1)%mp_test = 1
  allocate(mp_H_tmp(1)%Higgses(size(mp_H(index_i)%Higgses)+size(mp_H(index_j)%Higgses)))
  j=0
  do i=lbound(mp_H(index_i)%Higgses,dim=1),ubound(mp_H(index_i)%Higgses,dim=1)
   j=j+1
   mp_H_tmp(1)%Higgses(j)=mp_H(index_i)%Higgses(i)
  enddo
  do i=lbound(mp_H(index_j)%Higgses,dim=1),ubound(mp_H(index_j)%Higgses,dim=1)
   j=j+1
   mp_H_tmp(1)%Higgses(j)=mp_H(index_j)%Higgses(i)
  enddo
!-Now, fill with remaining, not combined Higgses
  j=1
  do i=lbound(mp_H,dim=1),ubound(mp_H,dim=1)
   if(mp_H(i)%mp_test.eq.1) then
    j=j+1
    allocate(mp_H_tmp(j)%Higgses(size(mp_H(i)%Higgses,dim=1)))
    mp_H_tmp(j)%Higgses = mp_H(i)%Higgses
    mp_H_tmp(j)%m = mp_H(i)%m
    mp_H_tmp(j)%dm = mp_H(i)%dm    
    mp_H_tmp(j)%mu = mp_H(i)%mu    
    mp_H_tmp(j)%mp_test = mp_H(i)%mp_test    
   endif
  enddo
 else
  finished_combination=1
 endif
   
end subroutine combine_Higgses 
!----------------------------------------------------------------------------
subroutine calculate_mpred_chisq(csq, N)
!----------------------------------------------------------------------------
 implicit none
 
 double precision, intent(out) :: csq
 integer, intent(out) :: N
 double precision, allocatable :: v(:), vmat(:,:), invcov(:,:), v2(:)
 double precision, allocatable :: csq_mu(:)
 integer :: i,ii,iii
 
 N=size(cov_mp,dim=1)
 allocate(v(N), vmat(N,1),invcov(N,N), v2(N), csq_mu(N))

 !-First construct the vector (mupred - muobs)_iii
 iii=0
 do i=lbound(analyses,dim=1),ubound(analyses,dim=1)
  do ii=lbound(analyses(i)%mpred%mp_Higgses,dim=1),&
&       ubound(analyses(i)%mpred%mp_Higgses,dim=1)
   iii=iii+1
    v(iii) = analyses(i)%mpred%mp_Higgses(ii)%mu-analyses(i)%mpred%mp_Higgses(ii)%mu_obs
   vmat(iii,1) = v(iii)
  enddo
 enddo
   
 call invmatrix(cov_mp,invcov)   
 call matmult(invcov,vmat,v2,N,1)
 
 iii=0
 do i=lbound(analyses,dim=1),ubound(analyses,dim=1)   
  do ii=lbound(analyses(i)%mpred%mp_Higgses,dim=1),&
&       ubound(analyses(i)%mpred%mp_Higgses,dim=1)
   iii=iii+1
   csq_mu(iii) = v(iii)*v2(iii)
   analyses(i)%mpred%mp_Higgses(ii)%chisq = csq_mu(iii)
  enddo
 enddo  

 csq=sum(csq_mu)
 deallocate(v,vmat,invcov,v2,csq_mu)


end subroutine calculate_mpred_chisq
!----------------------------------------------------------------------------
subroutine create_covariance_matrix_mp
!----------------------------------------------------------------------------
 implicit none

 integer :: N, i, ii, iii, j, jj, jjj
 double precision :: mumax, dmu0sq, dratesq
 character(LEN=50) :: title 
 
 if(.not.allocated(analyses)) then
  stop'Error in subroutine create_covariance_matrix_mp: analyses not allocated.' 
 endif
 if(allocated(cov_mp)) deallocate(cov_mp)

!! write(*,*) "Creating covariance matrix for mpred-method."

 iii=0
 do i=lbound(analyses,dim=1),ubound(analyses,dim=1)
  do ii=lbound(analyses(i)%mpred%mp_Higgses,dim=1),&
&       ubound(analyses(i)%mpred%mp_Higgses,dim=1)
!--Find out dimension (iii) of covariance matrix.
   iii=iii+1
   mumax=-1.0D10
!--Determine the dominating Higgs for every mpred-observable.
   do j=lbound(analyses(i)%mpred%mp_Higgses(ii)%Higgses,dim=1),&
&       ubound(analyses(i)%mpred%mp_Higgses(ii)%Higgses,dim=1)
    if(analyses(i)%mpred%mp_Higgses(ii)%Higgses(j)%mu.gt.mumax) then
     mumax=analyses(i)%mpred%mp_Higgses(ii)%Higgses(j)%mu
     analyses(i)%mpred%mp_Higgses(ii)%domH=analyses(i)%mpred%mp_Higgses(ii)%Higgses(j)%id
    endif     
   enddo
!--Subtract the correlated syst. uncertainties from (smeared) cyan band
   call correct_mu_uncertainty(analyses(i)%mpred%mp_Higgses(ii),analyses(i)%table)
!--(TS 08/05/2013) Calculate the channel weights in the model:
   call calculate_model_weights(analyses(i)%table, analyses(i)%mpred%mp_Higgses(ii))
!!   write(*,*) analyses(i)%mpred%mp_Higgses(ii)%channel_w_model
!!   write(*,*) analyses(i)%table%channel_w(:,analyses(i)%mpred%mp_Higgses(ii)%domH)
  enddo  
 enddo
 allocate(cov_mp(iii,iii))

 iii=0
 do i=lbound(analyses,dim=1),ubound(analyses,dim=1)
  do ii=lbound(analyses(i)%mpred%mp_Higgses,dim=1),&
&       ubound(analyses(i)%mpred%mp_Higgses,dim=1)
   iii=iii+1
   jjj=0
   do j=lbound(analyses,dim=1),ubound(analyses,dim=1)
    do jj=lbound(analyses(j)%mpred%mp_Higgses,dim=1),&
&         ubound(analyses(j)%mpred%mp_Higgses,dim=1)
     jjj=jjj+1
     cov_mp(iii,jjj)=0.0D0
     if(correlations_mu.or.(.not.correlations_mu.and.iii.eq.jjj)) then
      call get_rate_uncertainties_sq_mpred(dratesq, analyses(i)%table, &
&    analyses(i)%mpred%mp_Higgses(ii),analyses(j)%table,analyses(j)%mpred%mp_Higgses(jj))
!-----Treat luminosity uncertainty as correlated systematic error if same collaboration:
	  if(analyses(i)%table%collaboration.eq.analyses(i)%table%collaboration) then
       if(anticorrmu) then
 	   cov_mp(iii,jjj)=cov_mp(iii,jjj)+(analyses(i)%table%dlumi*analyses(j)%table%dlumi)&
& * (analyses(i)%mpred%mp_Higgses(ii)%mu*analyses(j)%mpred%mp_Higgses(jj)%mu)
       else
 	   cov_mp(iii,jjj)=cov_mp(iii,jjj)+(analyses(i)%table%dlumi*analyses(j)%table%dlumi)&       
& * abs(analyses(i)%mpred%mp_Higgses(ii)%mu*analyses(j)%mpred%mp_Higgses(jj)%mu)
        endif
	  endif
      if(anticorrmu) then
! (TS 28/10/2013: Change scaling of theo. err from observed to predicted mu:)
       cov_mp(iii,jjj)=cov_mp(iii,jjj)+(dratesq) * &
& (analyses(i)%mpred%mp_Higgses(ii)%total_mu*analyses(j)%mpred%mp_Higgses(jj)%total_mu)
      else
       cov_mp(iii,jjj)=cov_mp(iii,jjj)+(dratesq) * &
&  abs(analyses(i)%mpred%mp_Higgses(ii)%total_mu*analyses(j)%mpred%mp_Higgses(jj)%total_mu)
      endif
	 endif
!----Add the intrinsic (uncorrelated) uncertainty of this peak to the diagonal elements:	 
	 if(iii.eq.jjj) then
  	  call get_dmu0sq(dmu0sq,analyses(i)%mpred%mp_Higgses(ii))
	  cov_mp(iii,jjj) = cov_mp(iii,jjj) + dmu0sq
	 endif
    enddo 
   enddo  
  enddo
 enddo
 
end subroutine create_covariance_matrix_mp
!------------------------------------------------------------------------------------    
subroutine get_rate_uncertainties_sq_mpred(dratesq, table1, mp_H1, table2, mp_H2)
!------------------------------------------------------------------------------------    
 use usefulbits_HS, only : delta_rate
 
 type(mp_neutHiggs), intent(in) :: mp_H1, mp_H2
 type(mutable), intent(in) :: table1, table2
 double precision, intent(out) :: dratesq
 integer :: i,j,id1,p1,d1,id2,p2,d2
 double precision :: res
 res=0.0D0
!! write(*,*) 'table Nc = ',table1%Nc,table2%Nc
 do i=1,table1%Nc
  do j=1,table2%Nc
   id1 = table1%channel_id(i)
   p1 = int((id1-modulo(id1,10))/dble(10))
   d1 = modulo(id1,10)
   id2 = table2%channel_id(j)
   p2 = int((id2-modulo(id2,10))/dble(10))
   d2 = modulo(id2,10)
   if(p1.eq.p2) then
    res=res+delta_rate%dCS(p1)**2* &
&       mp_H1%channel_w_model(i)*mp_H2%channel_w_model(j)
   endif 
   if(d1.eq.d2) then
    res=res+delta_rate%dBR(d1)**2* &
&       mp_H1%channel_w_model(i)*mp_H2%channel_w_model(j)
   endif   
!   if(p1.eq.p2) then
!    res=res+delta_rate%dCS(p1)**2* &
!&           table1%channel_w(i,mp_H1%domH)*table2%channel_w(j,mp_H2%domH)
!   endif
!   if(d1.eq.d2) then
!    res=res+delta_rate%dBR(d1)**2* &
!&           table1%channel_w(i,mp_H1%domH)*table2%channel_w(j,mp_H2%domH)
!   endif
  enddo
 enddo
 
 dratesq=res
end subroutine get_rate_uncertainties_sq_mpred
!------------------------------------------------------------------------------------       
subroutine calculate_model_weights(table, mp_H)
 use usefulbits, only : small, div
 implicit none 
 
 type(mp_neutHiggs), intent(inout) :: mp_H
 type(mutable), intent(in) :: table

 integer :: k,l,Nc
 double precision :: normalization
 
!--In the (unphysical) case of negative channel rates, we have to take care that the
!  model channel weights are still positive and between 0 and 1:
   normalization = 0.0D0
   if(allocated(mp_H%channel_w_model)) deallocate(mp_H%channel_w_model)
   allocate(mp_H%channel_w_model(table%Nc))
   if(allocated(mp_H%channel_mu)) deallocate(mp_H%channel_mu)
   allocate(mp_H%channel_mu(table%Nc))


   mp_H%total_mu=0.0D0
   
   do k=lbound(mp_H%channel_w_model,dim=1),ubound(mp_H%channel_w_model,dim=1)
    mp_H%channel_mu(k)=0.0D0
    do l=lbound(mp_H%Higgses,dim=1),ubound(mp_H%Higgses,dim=1)
       mp_H%channel_mu(k)=mp_H%channel_mu(k) + &
&      table%channel_mu(k,mp_H%Higgses(l)%id)*table%channel_w(k,mp_H%Higgses(l)%id)
       mp_H%total_mu = mp_H%total_mu + &
&      table%channel_mu(k,mp_H%Higgses(l)%id)*table%channel_w(k,mp_H%Higgses(l)%id)
    enddo                     
   enddo


   normalization = 0.0D0
   do k=lbound(mp_H%channel_w_model,dim=1),ubound(mp_H%channel_w_model,dim=1)   
    mp_H%channel_w_model(k)=div(mp_H%channel_mu(k),mp_H%total_mu,0.0D0,1.0D9)
    normalization = normalization + abs(mp_H%channel_w_model(k))
   enddo

   if(abs(normalization).gt.small) then
    do k=lbound(mp_H%channel_w_model,dim=1),ubound(mp_H%channel_w_model,dim=1)
     mp_H%channel_w_model(k) = div(abs(mp_H%channel_w_model(k)),&
&                                  normalization,1.0D0/table%Nc,1.0D0)
    enddo
   else
!--If the predicted signal strength is zero, use SM weights   
    do k=lbound(mp_H%channel_w_model,dim=1),ubound(mp_H%channel_w_model,dim=1)
     mp_H%channel_w_model(k) = table%channel_w(k,mp_H%domH)
    enddo 
   endif 

   if(abs(sum(mp_H%channel_w_model)-1.0D0).gt.small) then
    write(*,*) "WARNING: Channel weights of the model are not correctly normalized:"
    write(*,*) mp_H%channel_w_model
   endif 


end subroutine calculate_model_weights


subroutine correct_mu_uncertainty(mp_H, table)
! The mu uncertainty as given in the mutable contains also systematic uncertainties for
! the luminosity and signal rate. These uncertainties are highly correlated to other
! analyses. Therefore, we subtract them here to obtain the intrinsic mu uncertainty of
! this peak. The correlated uncertainties enter later via the covariance matrix.
!------------------------------------------------------------------------------------    
 implicit none
 
 type(mp_neutHiggs), intent(inout) :: mp_H
 type(mutable), intent(in) :: table
 integer :: i
 double precision :: dcsq, dmuupsq, dmulowsq

 dcsq = 0
 do i=1, table%Nc
  dcsq = dcsq + table%channel_systSM(i,mp_H%domH)**2
 enddo

 dmuupsq = (mp_H%dmu_up_obs)**2-(table%dlumi*mp_H%mu_obs)**2-dcsq*mp_H%mu_obs**2
 dmulowsq = (mp_H%dmu_low_obs)**2-(table%dlumi*mp_H%mu_obs)**2-dcsq*mp_H%mu_obs**2

 mp_H%dmu_low0_obs = sqrt(abs((mp_H%dmu_low_obs)**2- &
&                             (table%dlumi*mp_H%mu_obs)**2-dcsq*mp_H%mu_obs**2))
 if(dmulowsq.lt.0.0D0)  mp_H%dmu_low0_obs = - mp_H%dmu_low0_obs
 mp_H%dmu_up0_obs = sqrt(abs((mp_H%dmu_up_obs)**2- &
&                            (table%dlumi*mp_H%mu_obs)**2-dcsq*mp_H%mu_obs**2))
 if(dmuupsq.lt.0.0D0)  mp_H%dmu_up0_obs = - mp_H%dmu_up0_obs

end subroutine correct_mu_uncertainty
!------------------------------------------------------------------------------------    
subroutine get_dmu0sq(dmu0sq,mp_H)
!------------------------------------------------------------------------------------    
 use usefulbits_hs, only : symmetricerrors
 implicit none
 
 double precision, intent(out) :: dmu0sq
 double precision :: pred_mu
 
 type(mp_neutHiggs), intent(in) :: mp_H 

 pred_mu = mp_H%mu

 if(.not.symmetricerrors) then 
  if(pred_mu.le.mp_H%mu_obs) then
   dmu0sq = mp_H%dmu_low0_obs**2
  else if(pred_mu.gt.mp_H%mu_obs) then
   dmu0sq = mp_H%dmu_up0_obs**2
  endif 
 else
  dmu0sq = (sqrt(mp_H%dmu_low0_obs**2)+sqrt(mp_H%dmu_up0_obs**2))**2/4.
 endif


end subroutine get_dmu0sq
!------------------------------------------------------------------------------------    
subroutine get_lower_or_upper_value_difference(val,pred,dataval)
!------------------------------------------------------------------------------------    
 implicit none
 
 double precision, intent(out) :: val
 double precision, intent(in) :: pred
 double precision, dimension(3) :: dataval
 
 if(pred.le.dataval(2)) then
  val = abs(dataval(2)-dataval(1))
 else if(pred.gt.dataval(2)) then
  val = abs(dataval(3)-dataval(2))
 endif 

end subroutine get_lower_or_upper_value_difference
!------------------------------------------------------------------------------------    
subroutine check_pdf
!------------------------------------------------------------------------------------    
 implicit none
 
 if(.not.((pdf.eq.1).or.(pdf.eq.2).or.(pdf.eq.3))) then
  write(*,*) 'WARNING: pdf not properly specified. Will be set to pdf=2 (gaussian-shape)'
  pdf=2
 endif
  
end subroutine check_pdf
!------------------------------------------------------------------------------------    
function Higgscomb_key(Higgses)
 implicit none
 type(neutHiggs), dimension(:), intent(in) :: Higgses
 integer :: Higgscomb_key
 integer :: i,keytmp,power

 keytmp=0
 power =0
 do i=lbound(Higgses,dim=1),ubound(Higgses,dim=1)
  keytmp = keytmp + Higgses(i)%id * 10**power
  power = power+1
 enddo
 Higgscomb_key = keytmp

end function Higgscomb_key
!------------------------------------------------------------------------------------
subroutine deallocate_mc_observables
 implicit none
 
 integer :: i,j
 do i=lbound(analyses,dim=1),ubound(analyses,dim=1)
  if(allocated(analyses(i)%mpred%mp_Higgses)) then
   do j=lbound(analyses(i)%mpred%mp_Higgses,dim=1),&
   &    ubound(analyses(i)%mpred%mp_Higgses,dim=1)
    if(allocated(analyses(i)%mpred%mp_Higgses(j)%Higgses)) then
     deallocate(analyses(i)%mpred%mp_Higgses(j)%Higgses)
    endif 
    if(allocated(analyses(i)%mpred%mp_Higgses(j)%channel_w_model)) then
     deallocate(analyses(i)%mpred%mp_Higgses(j)%channel_w_model)
    endif 
    if(allocated(analyses(i)%mpred%mp_Higgses(j)%channel_mu)) then
     deallocate(analyses(i)%mpred%mp_Higgses(j)%channel_mu)
    endif 
   enddo  
   deallocate(analyses(i)%mpred%mp_Higgses)
  endif
 enddo 

end subroutine deallocate_mc_observables
!------------------------------------------------------------------------------------
end module mc_chisq