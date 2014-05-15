!--------------------------------------------------------------------
! This file is part of HiggsSignals (TS 03/03/2013)
!--------------------------------------------------------------------
module all_chisq

 use numerics
 use combinatorics
 use usefulbits_hs
 implicit none
 
 contains
 
!--------------------------------------------------------------------
subroutine print_all_observables
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
  do j=lbound(analyses(i)%peaks,dim=1),ubound(analyses(i)%peaks,dim=1)
   write(*,'(A,1I3,A)') '#------------------------- peak observable ',j,&
&                       ' ----------------------------#'     
   write(*,'(A25,1F10.2)') 'peak_mass =',analyses(i)%peaks(j)%mpeak
   write(*,'(A25,1F10.4)') 'peak mu =',analyses(i)%peaks(j)%mu
   write(*,'(A25,2F10.4)') 'cyan band(low,high) =',analyses(i)%peaks(j)%dmulow,&
&    analyses(i)%peaks(j)%dmuup
   write(*,'(A25,4X,$)')   'Higgs combination ='
   do k=lbound(analyses(i)%peaks(j)%Higgs_comb,dim=1),&
&       ubound(analyses(i)%peaks(j)%Higgs_comb,dim=1)
    write(*,'(1I3,$)') analyses(i)%peaks(j)%Higgs_comb(k)
   enddo
   write(*,*) 
   write(*,'(A25,7X,1I3)') 'Dominant Higgs =',analyses(i)%peaks(j)%domH
   write(*,'(A25,1F10.4)'), 'Total pred. mu =',analyses(i)%peaks(j)%total_mu
   write(*,'(A25,1F15.9)'), 'Chisq for mu =',analyses(i)%peaks(j)%chisq_mu
   write(*,'(A25,1F15.9)'), 'Chisq for mh =',analyses(i)%peaks(j)%chisq_mh
   write(*,'(A25,1F15.9)'), 'Chisq (total) =',analyses(i)%peaks(j)%chisq_tot
   write(*,'(A25,1F15.9)'), 'Chisq (max) =',analyses(i)%peaks(j)%chisq_max
   write(*,*)'#------------------------ Channel information ----------------------------#'
   write(*,*)'  ID       prod.	decay 	       mu  	     weight  	 syst.err. '
   write(*,*)'#-------------------------------------------------------------------------#'   
   do k=1, analyses(i)%peaks(j)%Nc
    write(*,'(1I5,5X,2A,3F15.6)') analyses(i)%peaks(j)%channel_id(k),&
    &	analyses(i)%table%channel_description(k,:),&
    &	analyses(i)%peaks(j)%channel_mu(k),analyses(i)%peaks(j)%channel_w(k),&
    &	analyses(i)%peaks(j)%channel_syst(k)
   enddo 
  enddo    
  do j=lbound(analyses(i)%mpred%mp_Higgses,dim=1),&
&      ubound(analyses(i)%mpred%mp_Higgses,dim=1)
   write(*,*)'#--------------------- pred. mass centered observable --------------------#'   
   write(*,'(A25,4X,F10.2)') 'Higgs mass =', analyses(i)%mpred%mp_Higgses(j)%m
   write(*,'(A25,4X,F10.2)') 'Higgs mass uncertainty =',analyses(i)%mpred%mp_Higgses(j)%dm
   write(*,'(A25,4X,F10.6)') 'Pred. signal strength =',analyses(i)%mpred%mp_Higgses(j)%mu
   write(*,'(A25,4X,I10)') 'Combined Higgses =', &
&      size(analyses(i)%mpred%mp_Higgses(j)%Higgses,dim=1)
   do k=lbound(analyses(i)%mpred%mp_Higgses(j)%Higgses,dim=1),&
&       ubound(analyses(i)%mpred%mp_Higgses(j)%Higgses,dim=1)
   write(*,'(A5,I2,A18,4X,3F10.4)') 'Higgs ',k,' (m, dm, mu) =', &
&       analyses(i)%mpred%mp_Higgses(j)%Higgses(k)%m,  &
&       analyses(i)%mpred%mp_Higgses(j)%Higgses(k)%dm, &
&       analyses(i)%mpred%mp_Higgses(j)%Higgses(k)%mu
   enddo 
   write(*,'(A25,4X,F10.2)') 'obs. mass parameter =', &
&       analyses(i)%mpred%mp_Higgses(j)%m_obs
   write(*,'(A25,4X,F10.6)') 'obs. signal strength =', &
&       analyses(i)%mpred%mp_Higgses(j)%mu_obs
   write(*,'(A25,4X,2F10.6)') 'cyan band(low,high) =', &
&   analyses(i)%mpred%mp_Higgses(j)%dmu_low_obs,   &
&   analyses(i)%mpred%mp_Higgses(j)%dmu_up_obs
   write(*,'(A25,4X,2F10.6)') 'dmu0 (low,high) =', &
&   analyses(i)%mpred%mp_Higgses(j)%dmu_low0_obs,  &   
&   analyses(i)%mpred%mp_Higgses(j)%dmu_up0_obs
   write(*,'(A25,4X,F10.6)') 'Chisq =', analyses(i)%mpred%mp_Higgses(j)%chisq   
  enddo
 enddo 
 write(*,*) '#-------------------------------------------------------------------------#'   

end subroutine print_all_observables  
!------------------------------------------------------------------------------------    
subroutine calculate_total_chisq(csq_tot, csq_mu, csq_mh, csq_mpred, N, Nmu, Nmh, Nmpred)
! Still need to include assignment ranges here!
!------------------------------------------------------------------------------------    
 use pc_chisq, only : calculate_mh_chisq

 implicit none
 double precision, allocatable :: csq_mu_in(:),csq_mh_in(:,:),csq0(:),csq_mh_in_trunc(:,:)
 double precision, intent(out) :: csq_tot, csq_mu, csq_mh, csq_mpred
 integer, intent(out) :: N, Nmu, Nmh, Nmpred
 integer :: i, ii, iii, j, jjj, mpred_i

 call calculate_total_mu_chisq(csq0, Nmu, 1)
 call calculate_total_mu_chisq(csq_mu_in, Nmu, 0)	
 call calculate_mh_chisq(csq_mh_in, Nmh)

 allocate(csq_mh_in_trunc(size(csq_mh_in,dim=1),size(csq_mu_in,dim=1)))
 iii=0
 jjj=0
 do i=1,size(analyses)   
  do ii=1,size(analyses(i)%peaks,dim=1)+size(analyses(i)%mpred%mp_Higgses,dim=1) 
   if(ii.le.size(analyses(i)%peaks,dim=1)) then
    if(analyses(i)%peaks(ii)%NHiggs_comb.eq.0.and.&
&      size(analyses(i)%mpred%mp_Higgses).ne.0) then
     cycle
    endif 
   endif   
   jjj=jjj+1
   if(ii.le.size(analyses(i)%peaks,dim=1)) then
    iii=iii+1
    if(analyses(i)%peaks(ii)%NHiggs_comb.eq.0.and.&
&      size(analyses(i)%mpred%mp_Higgses).ne.0) then
	 N = N-1
     cycle
    endif 
    csq_mh_in_trunc(:,jjj) = csq_mh_in(:,iii)
   else
	csq_mh_in_trunc(:,jjj) = 0.0D0
   endif
  enddo
 enddo  

!-For debugging:
!! write(*,*) "csq_mh_in_trunc = "
 do i=lbound(csq_mh_in_trunc,dim=2),ubound(csq_mh_in_trunc,dim=2)
!!  write(*,*) csq_mh_in_trunc(:,i)
 enddo 
 
 csq_mu = 0.0D0 
 csq_mh = 0.0D0  
 csq_tot = 0.0D0
 csq_mpred = 0.0D0
 Nmpred = 0
 
 iii=0 
 do i=1,size(analyses)   
  do ii=1,size(analyses(i)%peaks,dim=1)+size(analyses(i)%mpred%mp_Higgses,dim=1)
   if(ii.le.size(analyses(i)%peaks,dim=1)) then
    if(analyses(i)%peaks(ii)%NHiggs_comb.eq.0.and.&
&      size(analyses(i)%mpred%mp_Higgses).ne.0) cycle
   endif 
   iii=iii+1  
!--Only allow positive chisq contributions if wanted   
   if(minimalchisq) then
    if(csq_mu_in(iii).lt.0.0D0) csq_mu_in(iii) = 0.0D0 
    do j=lbound(csq_mh_in_trunc, dim=1), ubound(csq_mh_in_trunc, dim=1)
     if(csq_mh_in_trunc(j,iii).lt.0.0D0) csq_mh_in_trunc(j,iii) = 0.0D0 
    enddo 
   endif 
!--Assign chisq_mu and chisq_mh such that the total chisq does not exceed the maximum
!--chisq csq0
   if(ii.le.size(analyses(i)%peaks,dim=1)) then
    if(maximalchisq) then	
     analyses(i)%peaks(ii)%chisq_mu = min(csq_mu_in(iii),csq0(iii))
     analyses(i)%peaks(ii)%chisq_mh = min(csq0(iii)-analyses(i)%peaks(ii)%chisq_mu,		&
&	 								     sum(csq_mh_in_trunc(:,iii)))
    elseif(analyses(i)%peaks(ii)%Higgs_assignment_forced.eq.0) then
     analyses(i)%peaks(ii)%chisq_mu = min(csq_mu_in(iii),csq0(iii))
     analyses(i)%peaks(ii)%chisq_mh = min(csq0(iii)-analyses(i)%peaks(ii)%chisq_mu,		&
&									    sum(csq_mh_in_trunc(:,iii)))
    else
     analyses(i)%peaks(ii)%chisq_mu = csq_mu_in(iii)
     analyses(i)%peaks(ii)%chisq_mh = sum(csq_mh_in_trunc(:,iii))
    endif
    analyses(i)%peaks(ii)%chisq_tot = analyses(i)%peaks(ii)%chisq_mu + 					&
& 									  analyses(i)%peaks(ii)%chisq_mh
    analyses(i)%peaks(ii)%chisq_max = csq0(iii)
    csq_mu = csq_mu + analyses(i)%peaks(ii)%chisq_mu
    csq_mh = csq_mh + analyses(i)%peaks(ii)%chisq_mh
    csq_tot = csq_tot + analyses(i)%peaks(ii)%chisq_tot
   else
    mpred_i = ii - size(analyses(i)%peaks,dim=1)
    analyses(i)%mpred%mp_Higgses(mpred_i)%chisq = csq_mu_in(iii)
    csq_mpred = csq_mpred + analyses(i)%mpred%mp_Higgses(mpred_i)%chisq
    csq_tot = csq_tot + analyses(i)%mpred%mp_Higgses(mpred_i)%chisq
    Nmpred = Nmpred + 1
   endif
  enddo 
 enddo 
 !Nmu is now only the number of signal strength peak observables:
 Nmu = Nmu - Nmpred
 !Calculate total ndf: 
 N = Nmu + Nmh + Nmpred

 deallocate(csq_mu_in, csq_mh_in, csq_mh_in_trunc)

end subroutine calculate_total_chisq
!----------------------------------------------------------------------------
subroutine calculate_total_mu_chisq(csq_mu, N, domax)
!----------------------------------------------------------------------------
 implicit none
 
 double precision, allocatable, intent(out) :: csq_mu(:)
 double precision :: csq
 integer, intent(out) :: N
 double precision, allocatable :: v(:), vmat(:,:), invcov(:,:), v2(:)
 integer :: i,ii,iii,mpred_i
 integer, intent(in) :: domax   ! if 1, then calculate maximal chisq

 if(allocated(cov_mu_tot)) deallocate(cov_mu_tot)
 call create_total_mu_covariance_matrix(domax)
 
 N=size(cov_mu_tot,dim=1)
 allocate(v(N), vmat(N,1),invcov(N,N), v2(N), csq_mu(N))

 !-First construct the vector (mupred - muobs)_iii
 iii=0
 do i=1,size(analyses)   
  do ii=1,size(analyses(i)%peaks,dim=1)+size(analyses(i)%mpred%mp_Higgses,dim=1)
   if(ii.le.size(analyses(i)%peaks,dim=1)) then
    if(analyses(i)%peaks(ii)%NHiggs_comb.eq.0.and.&
&      size(analyses(i)%mpred%mp_Higgses).ne.0) cycle
   endif 
   iii=iii+1
   if(ii.le.size(analyses(i)%peaks,dim=1)) then	
    if(domax.ge.1) then
     v(iii) = analyses(i)%peaks(ii)%mu
    else	
     v(iii) = analyses(i)%peaks(ii)%total_mu - analyses(i)%peaks(ii)%mu
    endif 
   else
    mpred_i = ii - size(analyses(i)%peaks,dim=1)
    v(iii) = analyses(i)%mpred%mp_Higgses(mpred_i)%mu - &
&            analyses(i)%mpred%mp_Higgses(mpred_i)%mu_obs
   endif 
   vmat(iii,1) = v(iii)
  enddo
 enddo
    
 call invmatrix(cov_mu_tot,invcov)   
 call matmult(invcov,vmat,v2,N,1)
 
 do i=1, N
  csq_mu(i) = v(i)*v2(i)
 enddo

 iii=0
 do i=1,size(analyses)   
  do ii=1,size(analyses(i)%peaks,dim=1)+size(analyses(i)%mpred%mp_Higgses,dim=1)
   if(ii.le.size(analyses(i)%peaks,dim=1)) then
    if(analyses(i)%peaks(ii)%NHiggs_comb.eq.0.and.&
&      size(analyses(i)%mpred%mp_Higgses).ne.0) cycle
   endif 
   iii=iii+1
   if(ii.le.size(analyses(i)%peaks,dim=1)) then
    if(domax.ge.1) then
     analyses(i)%peaks(ii)%chisq_max = csq_mu(iii)
    else 
     analyses(i)%peaks(ii)%chisq_mu = csq_mu(iii)
    endif 
   else
    mpred_i = ii - size(analyses(i)%peaks,dim=1)
    if(domax.eq.0) analyses(i)%mpred%mp_Higgses(mpred_i)%chisq = csq_mu(iii)    
   endif 
  enddo
 enddo

 csq=sum(csq_mu)
 deallocate(v,vmat,invcov,v2)

end subroutine calculate_total_mu_chisq
!----------------------------------------------------------------------------
subroutine create_total_mu_covariance_matrix(domax)
!------------------------------------------------------------------------------------
 use pc_chisq, only : get_dmu0sq_peak, get_rate_uncertainties_sq_peaks
 use mc_chisq, only : correct_mu_uncertainty, get_dmu0sq, get_rate_uncertainties_sq_mpred,&
&                     calculate_model_weights 
 implicit none
 integer, intent(in) :: domax   ! if 1, then use predicted mu == 0
 integer :: N, i, j, ii, jj, iii, jjj, mpred_i, mpred_j
 double precision :: mumax, dmu0sq, dratesq 
 character(LEN=50) :: title
 logical :: corr
 integer csqmax

!---TRY TO EVALUATE THE MAX CHISQ WITHOUT CORRELATIONS---!
 if(domax.eq.1) then
  csqmax = 1
  corr=correlations_mu
 elseif(domax.eq.2) then
  csqmax = 1
  corr=.False.
 else
  csqmax = 0
  corr=correlations_mu
 endif
!--------------------------------------------------------!  


!-First, find the dominant Higgs boson for each mpred observable and correct the signal strength
! uncertainty (i.e. subtract uncertainties which are added later as fully correlated uncertainties).
 iii=0
 do i=lbound(analyses,dim=1),ubound(analyses,dim=1)
  do ii=lbound(analyses(i)%mpred%mp_Higgses,dim=1),&
&       ubound(analyses(i)%mpred%mp_Higgses,dim=1)
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
  enddo
 enddo
 
!-Secondly, find the total number of (peak and mpred) observables, which have to be considered.
 N=0
 do i=lbound(analyses,dim=1),ubound(analyses,dim=1)
  do j=lbound(analyses(i)%peaks,dim=1),ubound(analyses(i)%peaks,dim=1)
   if(analyses(i)%peaks(j)%NHiggs_comb.eq.0.and.&
&     size(analyses(i)%mpred%mp_Higgses).eq.0)then
!---In this case, the observable cannot neither be covered by the peak nor the mpred method
!---, thus we want to take into account the penalty from the (peak) observable.
    N=N+1
   else if(analyses(i)%peaks(j)%NHiggs_comb.ne.0) then
!---Successfully assigned peak observable   
    N=N+1    
   endif 	
  enddo
!-Now, add also the number of mpred-observables for this analysis.
  N=N+size(analyses(i)%mpred%mp_Higgses)
 enddo
  
 allocate(cov_mu_tot(N,N))

 iii=0
 do i=1,size(analyses)   
  do ii=1,size(analyses(i)%peaks,dim=1)+size(analyses(i)%mpred%mp_Higgses,dim=1)
   if(ii.le.size(analyses(i)%peaks,dim=1)) then
    if(analyses(i)%peaks(ii)%NHiggs_comb.eq.0.and.&
&      size(analyses(i)%mpred%mp_Higgses).ne.0) cycle
   endif 
   iii=iii+1
   jjj=0
   do j=1, size(analyses)
    do jj=1,size(analyses(j)%peaks,dim=1)+size(analyses(j)%mpred%mp_Higgses,dim=1)
     if(jj.le.size(analyses(j)%peaks,dim=1)) then
      if(analyses(j)%peaks(jj)%NHiggs_comb.eq.0.and.&
&        size(analyses(j)%mpred%mp_Higgses).ne.0) cycle
     endif 
 	 jjj=jjj+1
	 cov_mu_tot(iii,jjj)=0.0D0
!----Treat the following cases: (1) peak x peak, (2) peak x mpred, (3) mpred x peak, (4) mpred x mpred
!-(1)------------------------------------------------------------------------------------
	 if(ii.le.size(analyses(i)%peaks,dim=1).and.jj.le.size(analyses(j)%peaks,dim=1)) then	 
  	  if(corr.or.(.not.corr.and.iii.eq.jjj)) then
      call get_rate_uncertainties_sq_peaks(dratesq,analyses(i)%peaks(ii),&
      &    analyses(j)%peaks(jj),analyses(i)%table%collider, analyses(j)%table%collider)  	  
!       call get_rate_uncertainties_sq_peaks(dratesq,analyses(i)%peaks(ii),&
!&                                           analyses(j)%peaks(jj))
       if(anticorrmu) then
        if(analyses(i)%table%collaboration.eq.analyses(j)%table%collaboration) then
 	     cov_mu_tot(iii,jjj)=cov_mu_tot(iii,jjj)+ &
&                           (analyses(i)%table%dlumi*analyses(j)%table%dlumi)* &
& 		  	  			    (analyses(i)%peaks(ii)%mu*analyses(j)%peaks(jj)%mu)
 	    endif
   	    cov_mu_tot(iii,jjj)=cov_mu_tot(iii,jjj)+(dratesq)* &
&                           (analyses(i)%peaks(ii)%total_mu*analyses(j)%peaks(jj)%total_mu)
       else
        if(analyses(i)%table%collaboration.eq.analyses(j)%table%collaboration) then
 	     cov_mu_tot(iii,jjj)=cov_mu_tot(iii,jjj)+ &
&                           (analyses(i)%table%dlumi*analyses(j)%table%dlumi)* &
& 		  	  			    abs(analyses(i)%peaks(ii)%mu*analyses(j)%peaks(jj)%mu)
 	    endif
   	    cov_mu_tot(iii,jjj)=cov_mu_tot(iii,jjj)+(dratesq)* &
&                          abs(analyses(i)%peaks(ii)%total_mu*analyses(j)%peaks(jj)%total_mu)
       endif        
 	  endif 
!----Add the intrinsic (uncorrelated) uncertainty of this peak to the diagonal elements:	 
	  if(iii.eq.jjj) then
  	   call get_dmu0sq_peak(dmu0sq,analyses(i)%peaks(ii),domax)
	   cov_mu_tot(iii,jjj) = cov_mu_tot(iii,jjj) + dmu0sq
	  endif
!-(2)------------------------------------------------------------------------------------
	 elseif(ii.le.size(analyses(i)%peaks,dim=1).and.&
&           jj.gt.size(analyses(j)%peaks,dim=1)) then
	  mpred_j = jj - size(analyses(j)%peaks,dim=1)
  	  if(corr) then
      call get_rate_uncertainties_sq_peak_mpred(dratesq, analyses(j)%table, &
&           analyses(j)%mpred%mp_Higgses(mpred_j), analyses(i)%peaks(ii))	   
       if(anticorrmu) then
   	    if(analyses(i)%table%collaboration.eq.analyses(j)%table%collaboration) then
 	     cov_mu_tot(iii,jjj)=cov_mu_tot(iii,jjj)+ &
&               (analyses(i)%table%dlumi*analyses(j)%table%dlumi)*  &
& 		  	  	(analyses(i)%peaks(ii)%mu*analyses(j)%mpred%mp_Higgses(mpred_j)%mu)
 	    endif
   	    cov_mu_tot(iii,jjj)=cov_mu_tot(iii,jjj)+(dratesq)*(analyses(i)%peaks(ii)%total_mu* &
&                          analyses(j)%mpred%mp_Higgses(mpred_j)%total_mu)
       else
   	    if(analyses(i)%table%collaboration.eq.analyses(j)%table%collaboration) then
 	     cov_mu_tot(iii,jjj)=cov_mu_tot(iii,jjj)+ &
&               (analyses(i)%table%dlumi*analyses(j)%table%dlumi)*  &
& 		  	  	abs(analyses(i)%peaks(ii)%mu*analyses(j)%mpred%mp_Higgses(mpred_j)%mu)
 	    endif
   	   cov_mu_tot(iii,jjj)=cov_mu_tot(iii,jjj)+(dratesq)*abs(analyses(i)%peaks(ii)%total_mu* &
&                          analyses(j)%mpred%mp_Higgses(mpred_j)%total_mu)
       endif
 	  endif 
!-(3)------------------------------------------------------------------------------------
	 elseif(ii.gt.size(analyses(i)%peaks,dim=1).and.&
&           jj.le.size(analyses(j)%peaks,dim=1)) then
	  mpred_i = ii - size(analyses(i)%peaks,dim=1)
  	  if(corr) then
      call get_rate_uncertainties_sq_peak_mpred(dratesq, analyses(i)%table, &
&           analyses(i)%mpred%mp_Higgses(mpred_i), analyses(j)%peaks(jj))	   
       if(anticorrmu) then
 	    if(analyses(i)%table%collaboration.eq.analyses(j)%table%collaboration) then
 	     cov_mu_tot(iii,jjj)=cov_mu_tot(iii,jjj)+ &
&           (analyses(i)%table%dlumi*analyses(j)%table%dlumi)* &
&           (analyses(j)%peaks(jj)%mu*analyses(i)%mpred%mp_Higgses(mpred_i)%mu)
 	    endif
   	    cov_mu_tot(iii,jjj)=cov_mu_tot(iii,jjj)+ &
&         (dratesq)*(analyses(j)%peaks(jj)%total_mu*analyses(i)%mpred%mp_Higgses(mpred_i)%total_mu)
       else
 	    if(analyses(i)%table%collaboration.eq.analyses(j)%table%collaboration) then
 	     cov_mu_tot(iii,jjj)=cov_mu_tot(iii,jjj)+ &
&           (analyses(i)%table%dlumi*analyses(j)%table%dlumi)* &
&           abs(analyses(j)%peaks(jj)%mu*analyses(i)%mpred%mp_Higgses(mpred_i)%mu)
 	    endif
   	    cov_mu_tot(iii,jjj)=cov_mu_tot(iii,jjj)+ &
&        (dratesq)*abs(analyses(j)%peaks(jj)%total_mu*analyses(i)%mpred%mp_Higgses(mpred_i)%total_mu)
       endif
 	  endif 
!-(4)------------------------------------------------------------------------------------
	 elseif(ii.gt.size(analyses(i)%peaks,dim=1).and.&
&           jj.gt.size(analyses(j)%peaks,dim=1)) then
	  mpred_i = ii - size(analyses(i)%peaks,dim=1)	 
	  mpred_j = jj - size(analyses(j)%peaks,dim=1)
     if(corr.or.(.not.corr.and.iii.eq.jjj)) then
      call get_rate_uncertainties_sq_mpred(dratesq, analyses(i)%table, &
&           analyses(i)%mpred%mp_Higgses(mpred_i),analyses(j)%table,   &
&           analyses(j)%mpred%mp_Higgses(mpred_j))
      if(anticorrmu) then
!-----Treat luminosity uncertainty as correlated systematic error if same collaboration:
 	   if(analyses(i)%table%collaboration.eq.analyses(i)%table%collaboration) then
 	    cov_mu_tot(iii,jjj)=cov_mu_tot(iii,jjj)+ &
&                          (analyses(i)%table%dlumi*analyses(j)%table%dlumi) * &
&                          (analyses(i)%mpred%mp_Higgses(mpred_i)%mu * &
&                           analyses(j)%mpred%mp_Higgses(mpred_j)%mu)
	   endif
   	   cov_mu_tot(iii,jjj)=cov_mu_tot(iii,jjj)+(dratesq)*								&
&  	   (analyses(i)%mpred%mp_Higgses(mpred_i)%total_mu*analyses(j)%mpred%mp_Higgses(mpred_j)%total_mu)	  
 	  else
!-----Treat luminosity uncertainty as correlated systematic error if same collaboration:
 	   if(analyses(i)%table%collaboration.eq.analyses(i)%table%collaboration) then
 	    cov_mu_tot(iii,jjj)=cov_mu_tot(iii,jjj)+ &
&                          (analyses(i)%table%dlumi*analyses(j)%table%dlumi) * &
&                          abs(analyses(i)%mpred%mp_Higgses(mpred_i)%mu * &
&                          analyses(j)%mpred%mp_Higgses(mpred_j)%mu)
	   endif
   	   cov_mu_tot(iii,jjj)=cov_mu_tot(iii,jjj)+(dratesq)*								&
&  	abs(analyses(i)%mpred%mp_Higgses(mpred_i)%total_mu*analyses(j)%mpred%mp_Higgses(mpred_j)%total_mu)
	  endif
	 endif 
!----Add the intrinsic (uncorrelated) uncertainty of this peak to the diagonal elements:	 
	 if(iii.eq.jjj) then
  	  call get_dmu0sq(dmu0sq,analyses(i)%mpred%mp_Higgses(mpred_i))
	  cov_mu_tot(iii,jjj) = cov_mu_tot(iii,jjj) + dmu0sq
	 endif
	 endif
  	enddo  	
   enddo 
  enddo 	 
 enddo 

!! title = "total (peak+mpred) covariance matrix for the signal strength"
!! call print_dble_matrix(cov_mu_tot,title)
 
end subroutine create_total_mu_covariance_matrix
!------------------------------------------------------------------------------------    
subroutine get_rate_uncertainties_sq_peak_mpred(dratesq, table, mp_H, peak)
!------------------------------------------------------------------------------------    
 implicit none
 type(mp_neutHiggs), intent(in) :: mp_H
 type(mutable), intent(in) :: table
 type(mupeak), intent(in) :: peak
 double precision, intent(out) :: dratesq
 integer :: i,j,id1,p1,d1,id2,p2,d2
 double precision :: res
 res=0.0D0
 do i=1,table%Nc
  do j=1,peak%Nc
   id1 = table%channel_id(i)
   p1 = int((id1-modulo(id1,10))/dble(10))
   d1 = modulo(id1,10)
   id2 = peak%channel_id(j)
   p2 = int((id2-modulo(id2,10))/dble(10))
   d2 = modulo(id2,10)

   if(p1.eq.p2) then
    res=res+delta_rate%dCS(p1)**2* &
&       mp_H%channel_w_model(i)**peak%channel_w_model(j)
   endif 
   if(d1.eq.d2) then
    res=res+delta_rate%dBR(d1)**2* &
&       mp_H%channel_w_model(i)**peak%channel_w_model(j)
   endif   
        
!   if(p1.eq.p2) res=res + &
!&                   delta_rate%dCS(p1)**2*table%channel_w(i,mp_H%domH)*peak%channel_w(j)
!   if(d1.eq.d2) res=res + &
!&                   delta_rate%dBR(d1)**2*table%channel_w(i,mp_H%domH)*peak%channel_w(j)
  enddo
 enddo
 
 dratesq=res
end subroutine get_rate_uncertainties_sq_peak_mpred
!------------------------------------------------------------------------------------    
end module all_chisq
!------------------------------------------------------------------------------------    