module expt_syst

! use usefulbits_hs
 implicit none

 integer,parameter :: Nprod = 5
! ggH, VBF, WH, ZH, ttH 
 integer,parameter :: Ndecay = 9
 integer,parameter :: Nsyst = 16
!  1: CMS H->gaga untagged 0-1 7 TeV event migration 
!  2: CMS H->gaga untagged 1-2 7 TeV event migration 
!  3: CMS H->gaga untagged 2-3 7 TeV event migration 
!  4: CMS H->gaga untagged 0-1 8 TeV event migration 
!  5: CMS H->gaga untagged 1-2 8 TeV event migration 
!  6: CMS H->gaga untagged 2-3 8 TeV event migration 
!  7: CMS H->gaga dijet 8 TeV event migration
!  8: Dijet tagging efficiency in dijet selection of CMS H->gaga analyses
!  9: ETmiss cut efficiency in ETmiss selection of CMS H->gaga analyses
! 10: ATLAS H->tautau ggH differential pT distribution and QCD scale
! 11: ATLAS H->tautau Top and Z->ll BG normalization (for hadlep and leplep channels)
! 12: ATLAS H->tautau hadronic tau identification and energy scale
! 13: ATLAS H->tautau JES eta calibration
! 14: ATLAS H->tautau Z->tautau normalization (for hadlep)
! 15: ATLAS H->tautau fake backgrounds (for leplep)
! 16: ATLAS H->tautau ditau(had) tagging efficiency
 double precision :: rel_corr_err(2,0:Nprod,0:Ndecay,Nsyst)
 
! scaletype determines whether the systematic uncertainty is scaled with the
! observed (typical for BG uncertainties) [0] or with the predicted mu (typical
! for signal uncertainties) [1].
 integer :: scaletype(Nsyst)
 
 contains
 
 !------------------------------------------------ 
 subroutine fill_scaletype
 !------------------------------------------------ 
 ! Set scaletypes (default is scaling with predicted) 
 scaletype(:)=1
 scaletype(11)=0
 scaletype(12)=0 
 scaletype(13)=0 
 scaletype(14)=0 
 scaletype(15)=0 
 scaletype(16)=0 
 
 end subroutine fill_scaletype
!------------------------------------------------ 
 subroutine fill_rel_corr_err(ID,N)
 !------------------------------------------------
  implicit none
   integer, intent(in) :: ID, N

  rel_corr_err(N,:,:,:) = 0.0D0
  
  select case(ID)

   case(13001107) ! untagged 0 7TeV
    rel_corr_err(N,:,1,1)=  +0.125D0
   case(13001108) ! untagged 1 7TeV
    rel_corr_err(N,:,1,1)=  -0.125D0
    rel_corr_err(N,:,1,2)=  +0.125D0
   case(13001109) ! untagged 2 7TeV
    rel_corr_err(N,:,1,2)=  -0.125D0
    rel_corr_err(N,:,1,3)=  +0.125D0    
   case(13001110) ! untagged 3 7TeV
    rel_corr_err(N,:,1,3)=  -0.125D0
   case(13001111) ! untagged 0 8TeV
    rel_corr_err(N,:,1,4)=  +0.125D0
   case(13001112) ! untagged 1 8TeV
    rel_corr_err(N,:,1,4)=  -0.125D0
    rel_corr_err(N,:,1,5)=  +0.125D0
   case(13001113) ! untagged 2 8TeV
    rel_corr_err(N,:,1,5)=  -0.125D0
    rel_corr_err(N,:,1,6)=  +0.125D0
   case(13001114) ! untagged 3 8TeV
    rel_corr_err(N,:,1,6)=  -0.125D0
   case(13001105) ! CMS H->gaga dijet loose tagged categories (8 TeV)
    rel_corr_err(N,:,1,7)= + 0.15D0
    rel_corr_err(N,1,1,8)= - 0.3D0
    rel_corr_err(N,2,1,8)= + 0.1D0
   case(13001106) ! CMS H->gaga dijet tight tagged categories (8 TeV)
    rel_corr_err(N,:,1,7)= - 0.15D0
    rel_corr_err(N,1,1,8)= - 0.3D0
    rel_corr_err(N,2,1,8)= + 0.15D0
   case(12015103) ! CMS H->gaga dijet tagged category (7 TeV)
    rel_corr_err(N,1,1,8)= - 0.3D0
    rel_corr_err(N,2,1,8)= + 0.1D0
   case(13001102) ! CMS H->gaga ETmiss tagged categories
    rel_corr_err(N,1,1,9)=   0.15D0 
    rel_corr_err(N,2,1,9)=   0.15D0
    rel_corr_err(N,3,1,9)=   0.04D0
    rel_corr_err(N,4,1,9)=   0.04D0
    rel_corr_err(N,5,1,9)=   0.04D0
   case(2012160101) ! ATL H->tautau leplep boosted category
    rel_corr_err(N,1,4,10)=   0.32D0
    rel_corr_err(N,:,4,11)=   0.15D0    
!    rel_corr_err(N,2,4,13)=  -0.12D0
    rel_corr_err(N,:,4,15)=   0.12D0
   case(2012160102) ! ATL H->tautau leplep VBF category 
    rel_corr_err(N,1,4,10)=   0.08D0
    rel_corr_err(N,:,4,11)=   0.15D0        
    rel_corr_err(N,:,4,13)=   0.12D0
!    rel_corr_err(N,2,4,13)=  -0.12D0
    rel_corr_err(N,:,4,15)=   0.12D0    
   case(2012160103) ! ATL H->tautau hadlep boosted category 
    rel_corr_err(N,1,4,10)=   0.32D0    
    rel_corr_err(N,:,4,11)=   0.15D0
    rel_corr_err(N,:,4,12)=   0.04D0
!    rel_corr_err(N,2,4,13)=  -0.12D0
    rel_corr_err(N,:,4,14)=   0.10D0    
   case(2012160104) ! ATL H->tautau hadlep VBF category    
    rel_corr_err(N,1,4,10)=   0.08D0
    rel_corr_err(N,:,4,11)=   0.15D0
    rel_corr_err(N,:,4,12)=   0.04D0!
    rel_corr_err(N,:,4,13)=   0.12D0
!    rel_corr_err(N,2,4,13)=  -0.12D0
    rel_corr_err(N,:,4,14)=   0.10D0
   case(2012160105) ! ATL H->tautau hadhad boosted category    
    rel_corr_err(N,1,4,10)=   0.22D0
    rel_corr_err(N,:,4,12)=   0.12D0
!    rel_corr_err(N,2,4,13)=  -0.12D0
    rel_corr_err(N,:,4,16)=   0.07D0
   case(2012160106) ! ATL H->tautau hadhad VBF category   
    rel_corr_err(N,1,4,10)=  0.05D0
    rel_corr_err(N,:,4,12)=  0.12D0
    rel_corr_err(N,:,4,13)=  0.12D0
!    rel_corr_err(N,2,4,13)= -0.12D0
    rel_corr_err(N,:,4,16)=  0.07D0        
   case default

  ! Good fit, but sloppy implementation:
!   case(13001107)
!    rel_corr_err(N,:,1,1)=   0.125D0
!   case(13001108)
!    rel_corr_err(N,:,1,1)=  -0.125D0
!   case(13001109)
!    rel_corr_err(N,:,1,1)=   0.125D0
!   case(13001110)
!    rel_corr_err(N,:,1,1)=  -0.125D0
!   case(13001111)
!    rel_corr_err(N,:,1,4)=   0.125D0
!   case(13001112)
!    rel_corr_err(N,:,1,4)=  -0.125D0
!   case(13001113)
!    rel_corr_err(N,:,1,4)=   0.125D0
!   case(13001114)
!    rel_corr_err(N,:,1,4)=  -0.125D0
!   case(13001105) ! CMS H->gaga dijet loose tagged categories (8 TeV)
!    rel_corr_err(N,1,1,2)= - 0.15D0
!    rel_corr_err(N,2,1,2)=   0.05D0
!   case(13001106) ! CMS H->gaga dijet tight tagged categories (8 TeV)
!    rel_corr_err(N,1,1,2)= + 0.3D0
!    rel_corr_err(N,2,1,2)= - 0.2D0
!   case(12015103) ! CMS H->gaga dijet tagged category (7 TeV)
!    rel_corr_err(N,1,1,2)= - 0.05D0
!    rel_corr_err(N,2,1,2)=   0.025D0
!   case(13001102) ! CMS H->gaga ETmiss tagged categories
!    rel_corr_err(N,1,1,3)=   0.15D0 
!    rel_corr_err(N,2,1,3)=   0.15D0
!    rel_corr_err(N,3,1,3)=   0.04D0
!    rel_corr_err(N,4,1,3)=   0.04D0
!    rel_corr_err(N,5,1,3)=   0.04D0
!   case default
  end select 
  
 end subroutine fill_rel_corr_err
!------------------------------------------------ 
 subroutine get_expt_syst_corr_for_peaks(value, peak1, mu1, peak2, mu2, model) 
!------------------------------------------------
  use usefulbits_hs, only : mupeak, print_dble_matrix
  implicit none
  
  type(mupeak), intent(in) :: peak1, peak2
  integer, intent(in) :: model
  double precision, intent(in) :: mu1, mu2 ! observed mu
  double precision, intent(out) :: value

  integer :: id1, p1, d1, id2, p2, d2  
  integer :: i,j,k
  
  call fill_rel_corr_err(peak1%id,1)
  call fill_rel_corr_err(peak2%id,2)  
  
!  if(peak1%id.eq.13001105.and.peak2%id.eq.13001106) then
!  write(*,*)'#-------------- ',peak1%id,' --------------#'
!  do k=1,Nprod
!   write(*,*) rel_corr_err(1,k,1,:)
!  enddo
!  write(*,*)'#-------------- ',peak2%id,' --------------#'
!  do k=1,Nprod
!   write(*,*) rel_corr_err(2,k,1,:)
!  enddo
!  endif
  
  value = 0.0D0  
  
  do i=lbound(peak1%channel_id,dim=1),ubound(peak1%channel_id,dim=1)
   do j=lbound(peak2%channel_id,dim=1),ubound(peak2%channel_id,dim=1)  
    id1 = peak1%channel_id(i)
    p1 = int((id1-modulo(id1,10))/dble(10))
    d1 = modulo(id1,10)
    id2 = peak2%channel_id(j)
    p2 = int((id2-modulo(id2,10))/dble(10))
    d2 = modulo(id2,10)

!  if(peak1%id.eq.13001105.and.peak2%id.eq.13001106) then
!   write(*,*) id1, p1, d1, id2, p2, d2
!   write(*,*) value
!  endif
    
    do k=1,Nsyst
     if(model.eq.1) then
      if(scaletype(k).eq.1) then
       value = value + &
&            peak1%channel_w_model(i)*rel_corr_err(1,p1,d1,k)*peak1%total_mu* &
&            peak2%channel_w_model(j)*rel_corr_err(2,p2,d2,k)*peak2%total_mu
      elseif(scaletype(k).eq.0) then
             value = value + &
&            peak1%channel_w(i)*rel_corr_err(1,p1,d1,k)*mu1* &
&            peak2%channel_w(j)*rel_corr_err(2,p2,d2,k)*mu2
      else
       write(*,*) "WARNING in get_expt_syst_corr_for peaks: Unknown scaletype of ",k
      endif 
     else
     value = value + &
&            peak1%channel_w(i)*rel_corr_err(1,p1,d1,k)*mu1* &
&            peak2%channel_w(j)*rel_corr_err(2,p2,d2,k)*mu2
     endif
    enddo
   enddo
  enddo 
  
!  if(abs(value).ge.0.0000001D0) then
!   write(*,*) "Non-zero correlated systematics:", peak1%id, peak2%id, mu1, mu2, value
!   write(*,*) "1st weights: ",peak1%channel_w_model
!   do k=1,Nprod
!    write(*,*) rel_corr_err(1,k,1,:)
!   enddo
!   write(*,*) "2nd weights: ",peak2%channel_w_model
!   do k=1,Nprod
!    write(*,*) rel_corr_err(2,k,1,:)
!   enddo
!
!  if(peak1%id.eq.13001105.and.peak2%id.eq.13001105) write(22,*) value
!  if(peak1%id.eq.13001106.and.peak2%id.eq.13001106) write(23,*) value, peak1%channel_w_model(1)*mu1
!  if(peak1%id.eq.12015103.and.peak2%id.eq.12015103) write(24,*) value, peak1%channel_w_model(1)*mu1
!  if(peak1%id.eq.13001105.and.peak2%id.eq.13001106) write(25,*) value
!  if(peak1%id.eq.12015103.and.peak2%id.eq.13001106) write(26,*) value


 end subroutine get_expt_syst_corr_for_peaks  
!------------------------------------------------  
end module  