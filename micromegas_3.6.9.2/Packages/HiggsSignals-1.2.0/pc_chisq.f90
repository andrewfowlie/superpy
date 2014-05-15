!--------------------------------------------------------------------
! This file is part of HiggsSignals (TS 03/03/2013)
!--------------------------------------------------------------------
module pc_chisq

 use numerics
 use combinatorics
 use usefulbits_hs
 implicit none

 integer :: i,j,k
 double precision,parameter :: pi=3.14159265358979323846264338328D0
 integer, allocatable :: peakindices_best(:,:)

contains

!------------------------------------------------------------------------------------    
subroutine set_toyvalues(ii, peaks)
! This subroutine sets the mass and mu measurements of the peak observable(s) of analyses ii
! to those values which are given by the user using assign_toyvalues_to_observables.
!------------------------------------------------------------------------------------   
 use usefulbits_hs, only: obs, mupeak
 
 integer, intent(in) :: ii
 type(mupeak),dimension(:), intent(inout) :: peaks
 integer :: i
 
 if(obs(ii)%table%npeaks.ne.size(peaks)) then
  stop 'Error in subroutine set_toyvalues: Number of peaks does not match!'
 endif
 
 do i=lbound(peaks,dim=1),ubound(peaks,dim=1)
  peaks(i)%mpeak = obs(ii)%table%Toys_mhobs(i)
  peaks(i)%mu = obs(ii)%table%Toys_muobs(i)
 enddo

end subroutine set_toyvalues
!------------------------------------------------------------------------------------   
!subroutine scale_uncertainties(ii, peaks)
!! Scales the uncertainty of the signal strength and mass measurement of the peak
!! observables of analysis ii by the scalefactors which have been set via the subroutine
!! assign_uncertainty_scalefactors_to_observables.
!!------------------------------------------------------------------------------------   
! use usefulbits_hs, only : obs, mupeak
!
! integer, intent(in) :: ii
! type(mupeak),dimension(:), intent(inout) :: peaks
! integer :: i
!
! if(obs(ii)%table%npeaks.ne.size(peaks)) then
!  stop 'Error in subroutine scale_uncertainties: Number of peaks does not match!'
! endif
! 
! do i=lbound(peaks,dim=1),ubound(peaks,dim=1)
!  peaks(i)%dmuup = obs(ii)%table%scale_mu(i)*peaks(i)%dmuup
!  peaks(i)%dmulow = obs(ii)%table%scale_mu(i)*peaks(i)%dmulow
!  peaks(i)%dm = obs(ii)%table%scale_mh(i)*peaks(i)%dm
! enddo
! 
!end subroutine scale_uncertainties 
!------------------------------------------------------------------------------------   
!subroutine restore_uncertainties(ii, peaks)
!! Restores the uncertainty of the signal strength and mass measurement of the peak
!! observables of analysis ii after scaling.
!!------------------------------------------------------------------------------------   
! use usefulbits, only : vsmall
! use usefulbits_hs, only : obs, mupeak
!
! integer, intent(in) :: ii
! type(mupeak),dimension(:), intent(inout) :: peaks
! integer :: i
!
! if(obs(ii)%table%npeaks.ne.size(peaks)) then
!  stop 'Error in subroutine restore_uncertainties: Number of peaks does not match!'
! endif
! 
! do i=lbound(peaks,dim=1),ubound(peaks,dim=1)
!  if(obs(ii)%table%scale_mu(i).ge.vsmall) then
!   peaks(i)%dmuup = peaks(i)%dmuup/obs(ii)%table%scale_mu(i)
!   peaks(i)%dmulow = peaks(i)%dmulow/obs(ii)%table%scale_mu(i)
!  else
!   write(*,*) 'WARNING: scale_mu is (too close to) zero!'
!  endif
!  if(obs(ii)%table%scale_mh(i).ge.vsmall) then
!   peaks(i)%dm = peaks(i)%dm/obs(ii)%table%scale_mh(i)
!  else
!   write(*,*) 'WARNING: scale_mh is (too close to) zero!'
!  endif   
! enddo
! 
!end subroutine restore_uncertainties 
!!------------------------------------------------------------------------------------   
subroutine assign_Higgs_to_peaks_with_correlations(iter)
! Do this only for pdf = 2
!
! NOTE: This is possibly still buggy. Only use it with iter=0.
! TODO: Have to extend this here for assignment-groups.
!------------------------------------------------------------------------------------   
 use usefulbits_HS, only : neutHiggses, nanalys
 use usefulbits, only : np, Hneut
 implicit none
 
 integer, intent(in) :: iter
 integer :: i, ii, iii, n, jjj
 character(LEN=100), allocatable :: assignmentgroups(:)
 integer, allocatable :: assignmentgroups_Higgs_comb(:,:)

!! allocate(assignmentgroups(nanalys))
!! allocate(assignmentgroups_Higgs_comb(nanalys,np(Hneut))) 

!! write(*,*) "Running assign_Higgs_to_peaks_with_correlations."

 if(iter.gt.0) write(*,*) "WARNING: Iterations in the Higgs-to-peaks assignment are ",&
&                         "still under development."

 jjj = 1
 do n=1, iter
  call create_covariance_matrices()
  iii=0
  do i=lbound(analyses,dim=1),ubound(analyses,dim=1)
   do ii=lbound(analyses(i)%peaks,dim=1),ubound(analyses(i)%peaks,dim=1)
    iii=iii+1
    analyses(i)%peaks(ii)%internalnumber = iii
   enddo 
   call assign_Higgs_to_peaks(analyses(i)%table, analyses(i)%peaks, n)
  enddo  
 enddo


end subroutine assign_Higgs_to_peaks_with_correlations
!------------------------------------------------------------------------------------    
subroutine assign_Higgs_to_peaks(table, peaks, iterstep)
! This subroutine assigns the best combination of Higgs bosons to each peak observable
! found in ONE mutable/analysis.
! It calculates for every possible assignment of the Higgs bosons to the peaks a
! chi-squared value. It takes care that each Higgs boson is used for at most one peak.
! The combination with the minimal chi-squared value is selected. The relevant information
! about the assigned Higgs bosons is then saved in each peak object.
!------------------------------------------------------------------------------------    
 use usefulbits, only : div, np, Hneut
 implicit none
 
 type(mutable), intent(in) :: table
 type(mupeak), dimension(:), intent(inout) :: peaks(:)
 integer, intent(in) :: iterstep

 integer, allocatable :: peakindices_best(:,:), domH(:), domH_tmp(:)
 integer, allocatable :: Higgs_to_be_assigned(:),Higgs_fulfills_req_assignment(:)
 integer :: Npeaks, nH, a, i, j, ii,jj, Hindex, NHiggs
 integer :: Higgs_assignment_forced,Higgs_assignment_forced_tmp
 double precision :: pccsq, pccsq_tmp, chisq_fp
!! integer, allocatable :: indices_best(:)
 double precision :: force_in_range

 chisq_fp = 10000000.0D0
 nH=np(Hneut)
 Npeaks=size(peaks)

 if(Npeaks.le.0) then 
  if(iterstep.eq.0) then
   write(*,'(A,1X,A,1X,A,1X,A,F4.2,A,I12,A)') ' No peaks defined for ', &
&   trim(adjustl(table%collaboration)),trim(adjustl(table%desc)),'(',table%energy, &
&   " TeV) search (analysis ID = ",table%id,")."
  endif 
 else
  !-Create indices matrix from combinatorics module containing all possible 
  !-peak-Higgs combinations:
  call create_peakindices(Npeaks,nH)
  allocate(peakindices_best(Npeaks,nH))
  allocate(domH(Npeaks),domH_tmp(Npeaks))
  allocate(Higgs_to_be_assigned(nH),Higgs_fulfills_req_assignment(nH))

! (TS 11/01/2013: Give some default values)
  Higgs_assignment_forced = 1
  call copy_matrices(peakindices(1,:,:),peakindices_best(:,:))
  do k=lbound(domH_tmp,dim=1),ubound(domH_tmp,dim=1)
   domH_tmp(k) = 0
  enddo  

 if(table%mhchisq.eq.1) then
  force_in_range = assignmentrange_massobs
 else 
  force_in_range = assignmentrange
 endif 
 if(pdf.eq.1) force_in_range = 1.0D0
  

! (TS 09/09/2012: Add the criterium that if the Higgs bosons masses are close enough to 
!  the peaks they have to be assigned. )
!-First, find out which Higgs bosons have to be assigned. For this, we loop over all peaks
!-and check which Higgses lie close enough (i.e. the mass difference is less than the
!-total (gaussian) mass uncertainty). In that case, we tag this Higgs to be assigned.
  do i=lbound(Higgs_to_be_assigned,dim=1), ubound(Higgs_to_be_assigned,dim=1)
   Higgs_to_be_assigned(i)=0
  enddo
  do i=lbound(peaks,dim=1),ubound(peaks,dim=1)
   do j=lbound(Higgs_to_be_assigned,dim=1), ubound(Higgs_to_be_assigned,dim=1)  
    if(abs(peaks(i)%Higgses(j)%m-peaks(i)%mpeak).le.									&
&	   force_in_range*sqrt(peaks(i)%Higgses(j)%dm**2 + table%deltam**2)) then
     Higgs_to_be_assigned(j)=1
    else
     !-If the chisq contribution from the Higgs masses is not used, we do NOT want
     ! to assign Higgs bosons which are far away from peak.
     if(table%mhchisq.ne.1) Higgs_to_be_assigned(j)=-1       
    endif 
   enddo
  enddo	

!-Loop over all possible combinations, calculate the chisq and find the best combination.
  do a=lbound(peakindices,dim=1),ubound(peakindices,dim=1)
   pccsq_tmp = 0
!--Loop over all peaks, i.e. rows of the combination matrices:
   do i=lbound(peakindices,dim=2),ubound(peakindices,dim=2)
!---Calculate the chi squared value for this combination assigned to peak i:
    call calc_pc_chisq(pccsq, peaks(i), table%mhchisq, peaks(i)%Higgses, domH(i),		&
&					   peakindices(a,i,:),iterstep)
!---Add the calculated value to the total chi-squared value for all peaks:	
    pccsq_tmp = pccsq_tmp + pccsq
   enddo

!--Determine the best Higgs-to-peaks assignment:	
   if(pccsq_tmp.lt.chisq_fp) then
     do i=lbound(Higgs_fulfills_req_assignment,dim=1), 									&
& 	      ubound(Higgs_fulfills_req_assignment,dim=1)
      Higgs_fulfills_req_assignment(i)=0
     enddo
     Higgs_assignment_forced_tmp = 0
     do j=lbound(Higgs_to_be_assigned,dim=1),ubound(Higgs_to_be_assigned,dim=1)
      if(Higgs_to_be_assigned(j).eq.1) then
       do ii=lbound(peakindices,dim=2),ubound(peakindices,dim=2)
        do jj=lbound(peakindices,dim=3),ubound(peakindices,dim=3)
         if(peakindices(a,ii,jj).eq.j) then
         Higgs_fulfills_req_assignment(j)=1
         Higgs_assignment_forced_tmp = 1
         endif
        enddo
       enddo
      else if(Higgs_to_be_assigned(j).eq.-1) then
       Higgs_fulfills_req_assignment(j)=1
       do ii=lbound(peakindices,dim=2),ubound(peakindices,dim=2)
        do jj=lbound(peakindices,dim=3),ubound(peakindices,dim=3)
         if(peakindices(a,ii,jj).eq.j) then
         Higgs_fulfills_req_assignment(j)=0
         endif
        enddo
       enddo
      else
       Higgs_fulfills_req_assignment(j)=1
      endif 
     enddo 
   	 if(sum(Higgs_fulfills_req_assignment).eq.nH) then
      chisq_fp = pccsq_tmp
      Higgs_assignment_forced = Higgs_assignment_forced_tmp
      call copy_matrices(peakindices(a,:,:),peakindices_best(:,:))
	  if(a.eq.1) then
	   do k=lbound(domH_tmp,dim=1),ubound(domH_tmp,dim=1)
        domH_tmp(k) = 0
	   enddo
	  else 
        domH_tmp(:) = domH(:)   
      endif  
     endif
   endif
  enddo
!--Save information in peak object
  do i=lbound(peakindices_best,dim=1),ubound(peakindices_best,dim=1)
!--Best Higgs Combination:
   peaks(i)%Higgs_comb(:) = peakindices_best(i,:)
!!   write(*,*) "hello: Higgs assignment of ID = ",peaks(i)%id, ": ", peaks(i)%Higgs_comb,&
!!   & "Higgs_assignment_forced: ",Higgs_assignment_forced
   peaks(i)%domH = domH_tmp(i)
   peaks(i)%Higgs_assignment_forced = Higgs_assignment_forced
   call evaluate_peak(peaks(i),table)
  enddo
  deallocate(peakindices, peakindices_best, domH, domH_tmp)
  deallocate(Higgs_to_be_assigned,Higgs_fulfills_req_assignment)  
  endif

end subroutine assign_Higgs_to_peaks
!------------------------------------------------------------------------------------    
subroutine evaluate_peak(peak, table)
! Evaluates the peak information for a given Higgs boson combination and dominant Higgs
! (both have to be assigned to the peak before calling this subroutine!)
!------------------------------------------------------------------------------------    
 use usefulbits, only : div, np, Hneut, small, vsmall
 implicit none
 type(mupeak) :: peak
 type(mutable) :: table
 integer :: j,k
 integer :: NHiggs,Hindex
 double precision :: normalization

 NHiggs=0
 do j=lbound(peak%Higgs_comb,dim=1),ubound(peak%Higgs_comb,dim=1)
  if(peak%Higgs_comb(j).ne.0) NHiggs=NHiggs+1
 enddo 
 peak%NHiggs_comb = NHiggs
!--Chose dominant Higgs in the best Higgs combination for weights and systematics 
!--------------
!  In rare cases there appears a segmentation fault, 
!  because peak%domH seems not to be initialized.
!  HERE: Check that it is in a reasonable range.
 if(peak%domH.ne.0.and.peak%domH.le.np(Hneut)) then
  peak%channel_w(:) = peak%channel_w_allH(:,peak%domH)
  peak%channel_w_corrected_eff(:) = peak%channel_w_corrected_eff_allH(:,peak%domH)  
  peak%channel_syst(:) = peak%channel_syst_allH(:,peak%domH)
  peak%channel_systSM(:) = peak%channel_systSM_allH(:,peak%domH)
 else
!-Part of fix: 
  peak%domH=0
!--------------  
  call get_weights_at_peak(peak, table)
!! If no Higgs is assigned we don't correct the efficiencies...  
  peak%channel_w_corrected_eff(:) = peak%channel_w(:)
!!  peak%channel_w(:) = peak%channel_w_allH(:,1)        
  peak%channel_syst(:) = peak%channel_syst_allH(:,1)     
  peak%channel_systSM(:) = peak%channel_systSM_allH(:,1)     													   
 endif
!--Subtract the correlated uncertainties from mu uncertainty:   
 call correct_mu_uncertainty_of_peak(peak, table)
!--Add the channel contributions of the Higgses and adjust them to the "averaged weights":
 peak%total_mu=0
!--n.b.: Have to set channel_mu to zero in case this subroutine is called several times.
 do k=lbound(peak%channel_mu,dim=1),ubound(peak%channel_mu,dim=1)
  peak%channel_mu(k)=0.0D0
  peak%channel_w_model(k)=0.0D0    
 enddo
!--Loop over Higgses in best combination
 do j=lbound(peak%Higgs_comb,dim=1),ubound(peak%Higgs_comb,dim=1)
  if(peak%Higgs_comb(j).ne.0) then
   Hindex=peak%Higgs_comb(j)
!----Loop over the channels and add rates
   do k=lbound(peak%channel_mu,dim=1),ubound(peak%channel_mu,dim=1)
    peak%channel_mu(k)=peak%channel_mu(k)+									&
&   peak%channel_mu_allH(k,Hindex)*peak%channel_w_corrected_eff_allH(k,Hindex)
    peak%total_mu = peak%total_mu + 											&
&	 peak%channel_mu_allH(k,Hindex)*peak%channel_w_corrected_eff_allH(k,Hindex)
   enddo     
  endif 
 enddo
   do k=lbound(peak%channel_mu,dim=1),ubound(peak%channel_mu,dim=1)
!---Calculate channel weights of the model, using the possibly combined rates.
!   (TS 20/04/2013)
    peak%channel_w_model(k)=div(peak%channel_mu(k),peak%total_mu,0.0D0,1.0D9)
!--Reweight the rates to obtain the channel_mu 
    peak%channel_mu(k)=div(peak%channel_mu(k),peak%channel_w_corrected_eff(k),0.0D0,1.0D9)
   enddo
!--If no Higgs boson has been assigned, use the SM weight at the peak position for
!  the channel weight of the model (TS26/04/2013):
   if(peak%domH.eq.0) then
    peak%channel_w_model = peak%channel_w
   endif
!--In the (unphysical) case of negative channel rates, we have to take care that the
!  model channel weights are still positive and between 0 and 1:
   normalization = 0.0D0
   do k=lbound(peak%channel_w_model,dim=1),ubound(peak%channel_w_model,dim=1)
    normalization = normalization + abs(peak%channel_w_model(k))
   enddo
   if(abs(normalization).gt.small) then
    do k=lbound(peak%channel_w_model,dim=1),ubound(peak%channel_w_model,dim=1)
     peak%channel_w_model(k) = div(abs(peak%channel_w_model(k)),&
&                                  normalization,1.0D0/peak%Nc,1.0D0)
    enddo
   else
!--If the predicted signal strength is zero, use SM weights   
    peak%channel_w_model = peak%channel_w 
   endif 
   if(abs(sum(peak%channel_w_model)-1.0D0).gt.small) then
    write(*,*) "WARNING: Channel weights of the model are not correctly normalized:"
    write(*,*) peak%channel_w_model
   endif 
end subroutine evaluate_peak
!------------------------------------------------------------------------------------    
subroutine calculate_total_pc_chisq(csq_tot, csq_mu, csq_mh, N, Nmu, Nmh)
! Calculates the total chi^2 value with the peak-centered chi^2 method. This 
! subroutine is called from evaluate_model in the HiggsSignals run.
!------------------------------------------------------------------------------------    
 implicit none
 double precision, allocatable :: csq_mu_in(:), csq_mh_in(:,:), csq0(:), csq0_nocorr(:)
 double precision, intent(out) :: csq_tot, csq_mu, csq_mh
 integer, intent(out) :: N, Nmu, Nmh
 integer :: i, ii, iii, j,jj, iter, Niterations
 logical :: iterate

! Need to do iterations of the following procedure. Due to modifications in the
! Higgs-to-peak assignments necessary for so-called assignment groups (i.e. a set of
! observables which should take over the assignment our their leading observable,
! which should be the one with a mass measurement), the covariance matrices of the Higgs
! mass part change, and therefore also the chi^2 contribution from the Higgs mass.
! Since a cutoff chi^2(tot) < chi^2(max) applies in specific cases, this cutoff has to be
! re-evaluated a few times to converge to the correct total chi^2 evaluation.

 iterate=.True.
 Niterations=0
 do while(iterate)

  call calculate_mu_chisq(csq0, Nmu, 2)
!  call calculate_mu_chisq(csq0_nocorr, Nmu, 2)      ! csq0 without correlations  
!  write(*,*) "Maximal chi2 values (i, with corr, without corr)"
!  do i=1,size(csq0)
!  write(*,*) i, csq0(i),  csq0_nocorr(i)
!  enddo
  call calculate_mu_chisq(csq_mu_in, Nmu, 0)	
  call calculate_mh_chisq(csq_mh_in, Nmh)
  
  iii=0 
  do i=1, size(analyses)
   do ii=lbound(analyses(i)%peaks,dim=1),ubound(analyses(i)%peaks,dim=1)
    iii=iii+1  
    if(analyses(i)%peaks(ii)%Higgs_assignment_forced.eq.0.and.(.not.maximalchisq)) then
     if((csq_mu_in(iii)+sum(csq_mh_in(:,iii))).gt.csq0(iii).and.&
&         analyses(i)%peaks(ii)%NHiggs_comb.gt.0.and.analyses(i)%table%mhchisq.eq.1) then
       analyses(i)%peaks(ii)%undo_assignment=1
!---Now, undo Higgs-to-peak-assignment for the whole group (if existent)	 
 	  if(len(trim(analyses(i)%peaks(ii)%assignmentgroup)).ne.0) then
 	   do j=1, size(analyses)
        do jj=lbound(analyses(j)%peaks,dim=1),ubound(analyses(j)%peaks,dim=1)
         if(analyses(i)%peaks(ii)%assignmentgroup.eq.&
&           analyses(j)%peaks(jj)%assignmentgroup) then
          analyses(j)%peaks(jj)%undo_assignment=1
         endif
        enddo
       enddo
      endif 
	 endif 
    endif
   enddo   
  enddo 

  call correcting_Higgs_to_peak_assignment(iterate)
  if(iterate) Niterations = Niterations + 1
  deallocate(csq_mu_in, csq_mh_in) 
 enddo
 
!! if(Niterations.gt.0) write(*,*) "Ran ",Niterations," iterations to determine Higgs to peak assignment."

 ! After these iterations, the code should know the correct assignments.
 call calculate_mh_chisq(csq_mh_in, Nmh)	! This will undo the assignments, if necessary
! call calculate_mu_chisq(csq0, Nmu, 1)
! call calculate_mu_chisq(csq0_nocorr, Nmu, 2)      ! csq0 without correlations
 call calculate_mu_chisq(csq_mu_in, Nmu, 0)	! Need to evaluate this again with new assignments
 
 !Calculate total ndf:
 N = Nmu + Nmh
 
 csq_mu = 0.0D0 
 csq_mh = 0.0D0  
 csq_tot = 0.0D0
 
 iii=0 
 do i=1, size(analyses)
  do ii=lbound(analyses(i)%peaks,dim=1),ubound(analyses(i)%peaks,dim=1)
   iii=iii+1  
!--Only allow positive chisq contributions if wanted   
   if(minimalchisq) then
    if(csq_mu_in(iii).lt.0.0D0) csq_mu_in(iii) = 0.0D0 
    do j=lbound(csq_mh_in, dim=1), ubound(csq_mh_in, dim=1)
     if(csq_mh_in(j,iii).lt.0.0D0) csq_mh_in(j,iii) = 0.0D0 
    enddo 
   endif 
!--Assign chisq_mu and chisq_mh such that the total chisq does not exceed the maximum
!--chisq csq0
   if(maximalchisq) then	
    analyses(i)%peaks(ii)%chisq_mu = min(csq_mu_in(iii),csq0(iii))
    analyses(i)%peaks(ii)%chisq_mh = min(csq0(iii)-analyses(i)%peaks(ii)%chisq_mu,		&
&									    sum(csq_mh_in(:,iii)))
!TESTING: comment out this -->
!   elseif(analyses(i)%peaks(ii)%Higgs_assignment_forced.eq.0) then
!    analyses(i)%peaks(ii)%chisq_mu = min(csq_mu_in(iii),csq0(iii))
!    analyses(i)%peaks(ii)%chisq_mh = min(csq0(iii)-analyses(i)%peaks(ii)%chisq_mu,		&
!&									    sum(csq_mh_in(:,iii)))
!!	write(*,*) csq0(iii),csq_mu_in(iii),sum(csq_mh_in(:,iii))
! NEW FOR CORRECTED ASSIGNMENTS:
!   elseif(analyses(i)%peaks(ii)%Higgs_assignment_forced.eq.1.and.&
!&         analyses(i)%peaks(ii)%domH.eq.0) then
    elseif(analyses(i)%peaks(ii)%domH.eq.0) then
    analyses(i)%peaks(ii)%chisq_mu = csq0(iii)
    analyses(i)%peaks(ii)%chisq_mh = 0.0D0
   else
!!    write(*,*) "HtP assignment is forced for analysis: ",iii   
    analyses(i)%peaks(ii)%chisq_mu = csq_mu_in(iii)
    analyses(i)%peaks(ii)%chisq_mh = sum(csq_mh_in(:,iii))
!!	write(*,*) csq0(iii),csq_mu_in(iii),sum(csq_mh_in(:,iii))    
   endif
   analyses(i)%peaks(ii)%chisq_tot = analyses(i)%peaks(ii)%chisq_mu + 					&
&									 analyses(i)%peaks(ii)%chisq_mh
   analyses(i)%peaks(ii)%chisq_max = csq0(iii)
   csq_mu = csq_mu + analyses(i)%peaks(ii)%chisq_mu
   csq_mh = csq_mh + analyses(i)%peaks(ii)%chisq_mh
   csq_tot = csq_tot + analyses(i)%peaks(ii)%chisq_tot
   
  enddo 
 enddo 

 deallocate(csq_mu_in, csq_mh_in)

end subroutine calculate_total_pc_chisq
!-----------------------------------------------------------------------------------
subroutine correcting_Higgs_to_peak_assignment(iterate)
!-----------------------------------------------------------------------------------
 use usefulbits, only : np, Hneut
 implicit none
 integer :: i, ii, k
 logical, intent(inout) :: iterate
 
 iterate=.False.
  do i=1, size(analyses)
   do ii=lbound(analyses(i)%peaks,dim=1),ubound(analyses(i)%peaks,dim=1)
    if(analyses(i)%peaks(ii)%undo_assignment.eq.1.and.analyses(i)%peaks(ii)%domH.ne.0)then
!!     write(*,*) "Correcting HtP."
     iterate=.True.
     do k=1,np(Hneut)
      analyses(i)%peaks(ii)%Higgs_comb(k)=0
     enddo
     analyses(i)%peaks(ii)%domH=0
     call evaluate_peak(analyses(i)%peaks(ii),analyses(i)%table)
!!     analyses(i)%peaks(ii)%undo_assignment=0
    endif
   enddo
  enddo 

end subroutine correcting_Higgs_to_peak_assignment
!------------------------------------------------------------------------------------
subroutine create_covariance_matrices()
!------------------------------------------------------------------------------------
 if(.not.allocated(analyses)) then
  stop 'Error in subroutine create_covariance_matrices: analyses is not allocated.'
 endif
 
 if(allocated(cov)) deallocate(cov)
 call create_covariance_matrix_mu(0)  
 if(pdf.eq.2) then
  if(allocated(cov_mhneut)) deallocate(cov_mhneut)
   call create_covariance_matrix_mhneut(0)
  if(allocated(cov_mhneut_max)) deallocate(cov_mhneut_max)
   call create_covariance_matrix_mhneut(1)
 endif

end subroutine create_covariance_matrices
!------------------------------------------------------------------------------------    
subroutine calculate_mu_chisq(csq_mu, N, domax)
!------------------------------------------------------------------------------------    
! use usefulbits_hs, only : peaklist, cov
 use numerics, only : invmatrix, matmult

 integer :: i, ii, iii
 double precision, allocatable :: v(:), vmat(:,:), invcov(:,:), v2(:)
 double precision, allocatable, intent(out) :: csq_mu(:)
 character(LEN=50) :: title = "covariance matrix for signal strength mu"
 integer, intent(out) :: N
 integer, intent(in) :: domax   ! if 1, then calculate maximal chisq
  
 if(allocated(cov)) deallocate(cov)
 call create_covariance_matrix_mu(domax) 
 
 if(allocated(mu_vector)) deallocate(mu_vector)

 N = size(cov,dim=1)
 allocate(v(N), vmat(N,1),invcov(N,N), v2(N), csq_mu(N))
 allocate(mu_vector(N))
 
!! write(*,*) ' domax = ',domax
 
 !-First construct the vector (mupred - muobs)_iii
 iii=0
 do i=1, size(analyses, dim=1)
  do ii=lbound(analyses(i)%peaks,dim=1),ubound(analyses(i)%peaks,dim=1)
   iii=iii+1
   if(domax.ge.1) then
    v(iii) = analyses(i)%peaks(ii)%mu
   else	
    v(iii) = analyses(i)%peaks(ii)%mu - analyses(i)%peaks(ii)%total_mu
   endif
   vmat(iii,1) = v(iii)
   
!!   write(*,*) 'v(',iii,') = ',v(iii)
  enddo
 enddo


! Copy vector into global module vector (for later access)
 mu_vector = v
    
 call invmatrix(cov,invcov)   
 call matmult(invcov,vmat,v2,N,1)
 
!! write(*,*) "Calculating mu chi^2. domax = ", domax
 
 do i=1, N
  csq_mu(i) = v(i)*v2(i)
!!  write(*,*) i, analyses(i)%peaks(1)%total_mu,csq_mu(i)
 enddo
!! write(*,*) "sum = ", sum(csq_mu(:))

 deallocate(v,vmat,invcov,v2)

end subroutine calculate_mu_chisq
!------------------------------------------------------------------------------------    
subroutine calculate_mh_chisq(csq_mh_out, ndf)
!------------------------------------------------------------------------------------    
! use usefulbits_hs, only : peaklist
 use usefulbits, only : np, Hneut
 use numerics, only : invmatrix, matmult

 integer, intent(out) :: ndf
 integer :: i, ii, iii, k, nH, N, Hindex
 double precision, allocatable :: v(:,:), invcov(:,:), v2(:), vmat(:,:)
 double precision, allocatable, intent(out) :: csq_mh_out(:,:)

 if((pdf.eq.1).or.(pdf.eq.3)) then
!-First, determine number of peaks: 
  N = 0
  do i=lbound(analyses,dim=1),ubound(analyses,dim=1)
   do j=lbound(analyses(i)%peaks,dim=1),ubound(analyses(i)%peaks,dim=1)
    N = N + 1
   enddo 
  enddo  
  nH = np(Hneut)
  allocate(csq_mh_out(nH,N))  
!-First, fill the chisq vector with zeros:
  do k=1,nH
   do i=1,N
    csq_mh_out(k,i) = 0.0D0
   enddo
  enddo   
  
  do k=1,nH
   iii=0
   do i=lbound(analyses,dim=1),ubound(analyses,dim=1)
    do ii=lbound(analyses(i)%peaks,dim=1),ubound(analyses(i)%peaks,dim=1)
     iii=iii+1
     Hindex = analyses(i)%peaks(ii)%Higgs_comb(k)
     if(Hindex.ne.0.and.analyses(i)%table%mhchisq.eq.1) then 
      csq_mh_out(Hindex,iii) = csq_mh(analyses(i)%peaks(ii)%Higgses(Hindex)%m, 			&
&        analyses(i)%peaks(ii)%mpeak,analyses(i)%peaks(ii)%Higgses(Hindex)%dm,			&
&		 analyses(i)%peaks(ii)%dm)
!!      write(*,*) 'chi2(mh',Hindex,')=',csq_mh_out(Hindex,iii)
!!      write(*,*) analyses(i)%peaks(ii)%dm, analyses(i)%peaks(ii)%Higgses(Hindex)%dm
	 endif
    enddo
   enddo 
  enddo 
  
 else if(pdf.eq.2) then
  if(allocated(cov_mhneut)) deallocate(cov_mhneut)
  call create_covariance_matrix_mhneut(0)

  nH = size(cov_mhneut,dim=1)
  N = size(cov_mhneut,dim=2)
  allocate(v(nH,N), v2(N), csq_mh_out(nH,N), vmat(N,1))
 
!-Construct the vector (mhpred - mhobs)_iii. Do this for every Higgs boson
!-of the model. If the Higgs boson (Hindex) is not assigned to the peak, we set
!-the entry of the vector to zero.
!-First, fill the vectors with zeros:
  do k=1,nH
   do i=1,N
    v(k,i) = 0.0D0
   enddo
  enddo  

  do k=1,nH
   iii=0
   do i=lbound(analyses,dim=1),ubound(analyses,dim=1)
    do ii=lbound(analyses(i)%peaks,dim=1),ubound(analyses(i)%peaks,dim=1)
     iii=iii+1
     Hindex = analyses(i)%peaks(ii)%Higgs_comb(k)
     if(Hindex.ne.0) then     
      v(Hindex,iii) = analyses(i)%peaks(ii)%Higgses(Hindex)%m -							&
&					  analyses(i)%peaks(ii)%mpeak
     endif  
    enddo
   enddo
  enddo 

  do k=1,nH 		!-n.b.: this loops now over Hindex
   call invmatrix(cov_mhneut(k,:,:),invcov) 
   vmat(:,1)=v(k,:)
!   call matmult(invcov,v(k,:),v2,N,1)
   call matmult(invcov,vmat,v2,N,1)   
   do i=1, N
    csq_mh_out(k,i) = v(k,i)*v2(i)
    if(csq_mh_out(k,i).ge.0.00001D0) then
!!     write(*,*) "hello: ",k,i,v(k,i),sign(1.0D0,v(k,i)), csq_mh_out(k,i)
    endif
   enddo   
  enddo  
  deallocate(v,v2,vmat,invcov)  
 endif

!--Determine number of observables ndf. This checks for each peak whether the mass
!--resolution is less than a huge number and thus the chisq evaluation from the mass
!--measurement is enabled.
  ndf=0
  if(pdf.eq.2.or.pdf.eq.3) then
   do i=lbound(analyses,dim=1),ubound(analyses,dim=1)
    do ii=lbound(analyses(i)%peaks,dim=1),ubound(analyses(i)%peaks,dim=1)
 	if(analyses(i)%table%mhchisq.eq.1) ndf=ndf+1
    enddo
   enddo 	
  endif

end subroutine calculate_mh_chisq
!------------------------------------------------------------------------------------    
subroutine create_covariance_matrix_mhneut(domax)
!------------------------------------------------------------------------------------    
 use usefulbits, only : np, Hneut
 implicit none

 double precision :: dratesq, dmu0sq, dmtemp
 character(LEN=50) :: title 	
 integer, intent(in) ::  domax   ! = 0 or 1. If set to 1, then fill all off-diagonal elements
 							! with squared theoretical mass uncertainty, corresponding to
 							! a complete Higgs-to-peaks assignment
 integer :: Npeaks,i,ii,iii,j,jj,jjj,k,Hindex_i, Hindex_j
  
!! write(*,*) "Create Higgs mass covariance matrix" 
  
 Npeaks = 0
 do i=lbound(analyses,dim=1),ubound(analyses,dim=1)
  do j=lbound(analyses(i)%peaks,dim=1),ubound(analyses(i)%peaks,dim=1)
   Npeaks = Npeaks + 1
  enddo 
 enddo
 
 if(domax.eq.0) then
  allocate(cov_mhneut(np(Hneut),Npeaks,Npeaks))
 elseif(domax.eq.1) then
  allocate(cov_mhneut_max(np(Hneut),Npeaks,Npeaks))
 else
  stop'ERROR in subroutine create_covariance_matrix_mhneut. Specify domax correctly!'
 endif
  
!-First, fill all elements of the covariance matrices with zero except the diagonal
!-elements which contain the experimental mass resolution squared
 do k=1,np(Hneut)
  iii=0
  do i=1, size(analyses)
   do ii=lbound(analyses(i)%peaks,dim=1),ubound(analyses(i)%peaks,dim=1)
    iii=iii+1
    jjj=0
    do j=1, size(analyses)
     do jj=lbound(analyses(j)%peaks,dim=1),ubound(analyses(j)%peaks,dim=1)
 	  jjj=jjj+1 	  
	  if(domax.eq.0) cov_mhneut(k,iii,jjj) = 0.0D0
	  if(domax.eq.1) cov_mhneut_max(k,iii,jjj) = 0.0D0
	  if(iii.eq.jjj) then
 !-----Deactivate mh chisq contributions for those analysis which are tagged mhchisq==0
 !-----by setting the uncertainty to a very large value. 
	   if(analyses(i)%table%mhchisq.eq.1) then
	    if(domax.eq.0) cov_mhneut(k,iii,jjj) = analyses(i)%peaks(ii)%dm**2
	    if(domax.eq.1) cov_mhneut_max(k,iii,jjj) = analyses(i)%peaks(ii)%dm**2
	   else
	    if(domax.eq.0) cov_mhneut(k,iii,jjj) = vlarge**2
	    if(domax.eq.1) cov_mhneut_max(k,iii,jjj) = vlarge**2	    
	   endif 	  
	  endif
     enddo
    enddo
   enddo 
  enddo
 enddo 
	     

!-Now, look for Higgs bosons shared by 2 peaks and fill the corresponding cov. matrix
!-element with the theoretical uncertainty squared of this Higgs boson.
 do k=1,np(Hneut)
  iii=0
  do i=1, size(analyses)
   do ii=lbound(analyses(i)%peaks,dim=1),ubound(analyses(i)%peaks,dim=1)
    iii=iii+1
    jjj=0
    do j=1, size(analyses)
     do jj=lbound(analyses(j)%peaks,dim=1),ubound(analyses(j)%peaks,dim=1)
 	  jjj=jjj+1 	  
      Hindex_i = analyses(i)%peaks(ii)%Higgs_comb(k)
      Hindex_j = analyses(j)%peaks(jj)%Higgs_comb(k)      
	  if((Hindex_i.eq.Hindex_j).and.(Hindex_i.ne.0).and.domax.eq.0) then
	   if(correlations_mh.or.(.not.correlations_mh.and.iii.eq.jjj)) then
        if(anticorrmh) then
         cov_mhneut(Hindex_i,iii,jjj)=cov_mhneut(Hindex_i,iii,jjj) + 					&
& sign(1.0D0,(analyses(i)%peaks(ii)%Higgses(Hindex_i)%m-analyses(i)%peaks(ii)%mpeak))*  &
& sign(1.0D0,(analyses(j)%peaks(jj)%Higgses(Hindex_j)%m-analyses(j)%peaks(jj)%mpeak))*  &
& analyses(i)%peaks(ii)%Higgses(Hindex_i)%dm**2
         else
          cov_mhneut(Hindex_i,iii,jjj)=cov_mhneut(Hindex_i,iii,jjj) + 					&
& analyses(i)%peaks(ii)%Higgses(Hindex_i)%dm**2 
         endif
   	   endif
	  elseif(domax.eq.1) then
	   if(correlations_mh.or.(.not.correlations_mh.and.iii.eq.jjj)) then
!	    
!NEW TEST (sign dependence):
        if(anticorrmh) then
         cov_mhneut_max(k,iii,jjj)= cov_mhneut_max(k,iii,jjj) + 							&
& sign(1.0D0,(analyses(i)%peaks(ii)%Higgses(k)%m-analyses(i)%peaks(ii)%mpeak))*  &
& sign(1.0D0,(analyses(j)%peaks(jj)%Higgses(k)%m-analyses(j)%peaks(jj)%mpeak))*  &
& analyses(i)%peaks(ii)%Higgses(k)%dm**2
        else
        cov_mhneut_max(k,iii,jjj)= cov_mhneut_max(k,iii,jjj) + 							&
& analyses(i)%peaks(ii)%Higgses(k)%dm**2
        endif
	   endif
  	  endif  	  
     enddo
    enddo 
   enddo 	 
  enddo 
 enddo
!! title = "mass covariance matrix for maximal assignment, first Higgs boson"
!! if(domax.eq.1) call print_dble_matrix(cov_mhneut_max(1,:,:), title)
!! title = "covariance matrix for signal strength mh2"
!! call print_dble_matrix(cov_mhneut(2,:,:), title)
!! title = "covariance matrix for signal strength mh3" 
!! call print_dble_matrix(cov_mhneut(3,:,:), title) 
end subroutine create_covariance_matrix_mhneut
!------------------------------------------------------------------------------------ 
subroutine correct_mu_uncertainty_of_peak(peak, table)
! The mu uncertainty as given in the mutable contains also systematic uncertainties for
! the luminosity and signal rate. These uncertainties are highly correlated to other
! analyses. Therefore, we subtract them here to obtain the intrinsic mu uncertainty of
! this peak. The correlated uncertainties enter later via the covariance matrix.
!------------------------------------------------------------------------------------    
 use usefulbits_hs, only : usescalefactor, output_level,symmetricerrors,&
& absolute_errors, THU_included!, withcorrexpsyst
 implicit none

 type(mupeak), intent(inout) :: peak
 type(mutable), intent(in) :: table
! double precision, intent(in) :: dlumi  	! Uncertainty of Luminosity (relative)
 integer :: j
 double precision :: dcsq, dmulow0sq, dmuup0sq, allsystsq, dmuaverage, dmuaverage0sq, mu

 dcsq = 0
 do j=1, peak%Nc
  dcsq = dcsq + peak%channel_systSM(j)**2
 enddo

 if(absolute_errors) then
  mu=peak%mu_original
 else
  mu=peak%mu
 endif

! if(.not.symmetricerrors) then
 if(THU_included) then
!  if(withcorrexpsyst) then
!   allsystsq = 0.0D0
!   do j=lbound(table%correxpsyst,dim=1),ubound(table%correxpsyst,dim=1)
!    allsystsq = allsystsq + table%correxpsyst(j)**2
!   enddo 
!   dmulow0sq = (peak%dmulow)**2-(allsystsq + dcsq)*mu**2
!   dmuup0sq = (peak%dmuup)**2-(allsystsq + dcsq)*mu**2
!  else 
   dmulow0sq = (peak%dmulow)**2-(table%dlumi*mu)**2-dcsq*mu**2
   dmuup0sq = (peak%dmuup)**2-(table%dlumi*mu)**2-dcsq*mu**2
!  endif
 else
!- In the case of future projections we usually use measurements
!- without theoretical uncertainty included:
  dmulow0sq = (peak%dmulow)**2
  dmuup0sq = (peak%dmuup)**2
 endif
  
  if(usescalefactor) then
   dmulow0sq = peak%scale_mu**2*dmulow0sq
   dmuup0sq =  peak%scale_mu**2*dmuup0sq
  endif
  if(.not.symmetricerrors) then
  peak%dmulow0sq = dmulow0sq
  peak%dmuup0sq = dmuup0sq
  else
   if(peak%dmulow0sq.lt.0.0d0.or.peak%dmuup0sq.lt.0.0d0) then
    write(*,*) "WARNING: squared intrinsic mu uncertainty is negative!"
   endif 
   peak%dmulow0sq = (sqrt(abs(dmulow0sq))+sqrt(abs(dmuup0sq)))**2/4.0d0
   peak%dmuup0sq = (sqrt(abs(dmulow0sq))+sqrt(abs(dmuup0sq)))**2/4.0d0
  endif 
! else
!  dmuaverage = (peak%dmulow+peak%dmuup)/2.
!  if(withcorrexpsyst) then
!  allsystsq = 0.0D0
!   do j=lbound(table%correxpsyst,dim=1),ubound(table%correxpsyst,dim=1)
!    allsystsq = allsystsq + table%correxpsyst(j)**2
!   enddo 
!   dmuaverage0sq = (dmuaverage)**2-(allsystsq + dcsq)*peak%mu**2
!  else 
!   dmuaverage0sq = (dmuaverage)**2-(table%dlumi*peak%mu)**2-dcsq*peak%mu**2
!  endif
!  if(usescalefactor) then
!   dmuaverage0sq = peak%scale_mu**2*dmuaverage0sq
!  endif
!  peak%dmulow0sq = dmuaverage0sq
!  peak%dmuup0sq = dmuaverage0sq
! endif
     


! if(output_level.ne.0) then
  if(peak%dmulow0sq.lt.0.0D0) then
   write(*,*) "WARNING: Negative intrinsic (lower) error squared for observable ID = ",&
   			 peak%ID,", with value ",peak%dmulow0sq
  endif 
  if(peak%dmuup0sq.lt.0.0D0) then
   write(*,*) "WARNING: Negative intrinsic (upper) error squared for observable ID = ",&
   			 peak%ID,", with value ",peak%dmuup0sq
  endif 
! endif 

 
!! write(*,*) "Original / Corrected uncertainties (low): ",peak%dmulow, sqrt(dmulow0sq)
!! write(*,*) "Original / Corrected uncertainties (up): ",peak%dmuup, sqrt(dmuup0sq)
 
 
end subroutine correct_mu_uncertainty_of_peak
!------------------------------------------------------------------------------------    
subroutine create_covariance_matrix_mu(domax)
!------------------------------------------------------------------------------------    
 use usefulbits_hs, only : withcorrexpsyst, absolute_errors
 use expt_syst, only : get_expt_syst_corr_for_peaks
 implicit none

 integer, intent(in) :: domax   ! if 1, then use predicted mu == 0
 double precision :: dratesq, dmu0sq, s1, s2, mu_iii, mu_jjj, corrsyst, corrsystSM
 integer :: Npeaks
 integer :: i,ii,iii,j,jj,jjj
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
  
 Npeaks = 0
 do i=lbound(analyses,dim=1),ubound(analyses,dim=1)
   Npeaks = Npeaks + analyses(i)%Npeaks
 enddo
   
 allocate(cov(Npeaks,Npeaks))
 
 iii=0
 !--Loop twice over all peaks to construct matrix
 do i=1, size(analyses)
  do ii=lbound(analyses(i)%peaks,dim=1),ubound(analyses(i)%peaks,dim=1)
   iii=iii+1
   if(absolute_errors) then
    mu_iii=analyses(i)%peaks(ii)%mu_original
   else
    mu_iii=analyses(i)%peaks(ii)%mu
   endif   
   jjj=0
   do j=1, size(analyses)
    do jj=lbound(analyses(j)%peaks,dim=1),ubound(analyses(j)%peaks,dim=1)
 	 jjj=jjj+1
     if(absolute_errors) then
      mu_jjj=analyses(j)%peaks(jj)%mu_original
     else
      mu_jjj=analyses(j)%peaks(jj)%mu
     endif

 	 if(usescalefactor) then
 	  s1 = analyses(i)%peaks(ii)%scale_mu
   	  s2 = analyses(j)%peaks(jj)%scale_mu
   	 else
   	  s1 = 1.0D0
   	  s2 = 1.0D0
   	 endif  
	 cov(iii,jjj)=0.0D0
 	 if(corr.or.(.not.corr.and.iii.eq.jjj)) then
      call get_rate_uncertainties_sq_peaks(dratesq,analyses(i)%peaks(ii),&
      &    analyses(j)%peaks(jj),analyses(i)%table%collider, analyses(j)%table%collider)
      if(anticorrmu) then
!-----Treat luminosity uncertainty as correlated systematic error if same collaboration:
 	   if(analyses(i)%table%collaboration.eq.analyses(j)%table%collaboration) then
!!        if(withcorrexpsyst) then !!THIS IS IN TESTING PHASE!!
!!         do k=lbound(analyses(i)%table%correxpsyst,dim=1),&
!!&            ubound(analyses(i)%table%correxpsyst,dim=1)
!!         cov(iii,jjj)=cov(iii,jjj) + ( s1*analyses(i)%table%correxpsyst(k)  &
!!&                                    *s2*analyses(j)%table%correxpsyst(k) ) &
!!& 		 	  					  * (mu_iii*mu_jjj)
!!         enddo
!!        else
  	     cov(iii,jjj)=cov(iii,jjj)+(s1*analyses(i)%table%dlumi*s2*analyses(j)%table%dlumi)&
& 		 	  					*(mu_iii*mu_jjj)
!!        endif
	   endif
	   corrsyst=0.0D0

!------Include correlated experimental systematic uncertainties (TS 2013-11-21)
       if(withcorrexpsyst) then       
        if(iii.eq.jjj) then
         call get_expt_syst_corr_for_peaks(corrsystSM, analyses(i)%peaks(ii),&
&             mu_iii, analyses(j)%peaks(jj),mu_jjj, 0) 
         cov(iii,jjj)=cov(iii,jjj) - corrsystSM
        endif
!! TESTING: 21/01/2014
        call get_expt_syst_corr_for_peaks(corrsyst, analyses(i)%peaks(ii),&
&            mu_iii, analyses(j)%peaks(jj),mu_jjj, 1) 
!!
!!        call get_expt_syst_corr_for_peaks(corrsyst, analyses(i)%peaks(ii),&
!!&            analyses(i)%peaks(ii)%total_mu, analyses(j)%peaks(jj),analyses(j)%peaks(jj)%total_mu, 1) 
        cov(iii,jjj)=cov(iii,jjj) + corrsyst
       endif 
!------
       if(useSMweights) then
        cov(iii,jjj)=cov(iii,jjj)+ &
&       (dratesq)*(analyses(i)%peaks(ii)%mu*analyses(j)%peaks(jj)%mu)
       else 
   	    cov(iii,jjj)=cov(iii,jjj)+ &
&       (dratesq)*(analyses(i)%peaks(ii)%total_mu*analyses(j)%peaks(jj)%total_mu)
       endif
      else
!-----Treat luminosity uncertainty as correlated systematic error if same collaboration:
 	   if(analyses(i)%table%collaboration.eq.analyses(j)%table%collaboration) then
!!        if(withcorrexpsyst) then !!THIS IS IN TESTING PHASE!!
!!         do k=lbound(analyses(i)%table%correxpsyst,dim=1),&
!!&            ubound(analyses(i)%table%correxpsyst,dim=1)
!!         cov(iii,jjj)=cov(iii,jjj) + ( s1*analyses(i)%table%correxpsyst(k)  &
!!&                                    *s2*analyses(j)%table%correxpsyst(k) )&
!!& 		 	  					  * abs(mu_iii*mu_jjj)
!!         enddo
!!        else
  	     cov(iii,jjj)=cov(iii,jjj)+(s1*analyses(i)%table%dlumi*s2*analyses(j)%table%dlumi)&
& 		 	  					  *abs(mu_iii*mu_jjj)
!!        endif
	   endif
       if(useSMweights) then
        cov(iii,jjj)=cov(iii,jjj)+ &
&       (dratesq)*(analyses(i)%peaks(ii)%mu*analyses(j)%peaks(jj)%mu)
       else 
   	    cov(iii,jjj)=cov(iii,jjj)+ &
&       (dratesq)*(analyses(i)%peaks(ii)%total_mu*analyses(j)%peaks(jj)%total_mu)
       endif
      endif        
 	 endif 
!----Add the intrinsic (uncorrelated) uncertainty of this peak to the diagonal elements:	 
	 if(iii.eq.jjj) then
  	  call get_dmu0sq_peak(dmu0sq,analyses(i)%peaks(ii),csqmax)
	  cov(iii,jjj) = cov(iii,jjj) + dmu0sq
	 endif
  	enddo  	
   enddo 
  enddo 	 
 enddo 
 
!! write(*,*) "signal strength covariance matrix:"
!! do iii=lbound(cov,dim=1),ubound(cov,dim=1)
!!  write(*,*) "(",iii,",",iii,") = ",cov(iii,iii)
!! enddo
end subroutine create_covariance_matrix_mu
!------------------------------------------------------------------------------------    
subroutine get_rate_uncertainties_sq_peaks_old(dratesq, peak1, peak2)
! Returns the sum of the squared theoretical uncertainties of the common 
! production and decay rates of the two peak observables.
!------------------------------------------------------------------------------------    
 use usefulbits_HS, only : mupeak, delta_rate, useSMweights
 
 type(mupeak), intent(in) :: peak1, peak2
 double precision, intent(out) :: dratesq
 integer :: i,j,id1,p1,d1,id2,p2,d2
 double precision :: res
 res=0.0D0
 do i=1,peak1%Nc
  do j=1,peak2%Nc
   id1 = peak1%channel_id(i)
   p1 = int((id1-modulo(id1,10))/dble(10))
   d1 = modulo(id1,10)
   id2 = peak2%channel_id(j)
   p2 = int((id2-modulo(id2,10))/dble(10))
   d2 = modulo(id2,10)
!   if(p1.eq.p2) res=res+delta_rate%dCS(p1)**2*peak1%channel_w(i)*peak2%channel_w(j)
!   if(d1.eq.d2) res=res+delta_rate%dBR(d1)**2*peak1%channel_w(i)*peak2%channel_w(j)
   if(p1.eq.p2) then
    if(useSMweights) then
     res=res+delta_rate%dCS(p1)**2*peak1%channel_w(i)*peak2%channel_w(j)
    else
     res=res+delta_rate%dCS(p1)**2*peak1%channel_w_model(i)*peak2%channel_w_model(j)
    endif 
   endif 
   if(d1.eq.d2) then
    if(useSMweights) then
     res=res+delta_rate%dBR(d1)**2*peak1%channel_w(i)*peak2%channel_w(j)
    else 
     res=res+delta_rate%dBR(d1)**2*peak1%channel_w_model(i)*peak2%channel_w_model(j)
    endif 
   endif   
  enddo
 enddo
  
 dratesq=res
end subroutine get_rate_uncertainties_sq_peaks_old
!------------------------------------------------------------------------------------    
subroutine get_rate_uncertainties_sq_peaks(dratesq, peak1, peak2, collider1, collider2)
! Returns the sum of the squared theoretical uncertainties of the common 
! production and decay rates of the two peak observables.
!------------------------------------------------------------------------------------    
 use usefulbits_HS, only : mupeak, delta_rate, useSMweights
 
 type(mupeak), intent(in) :: peak1, peak2
 character(LEN=*), intent(in) :: collider1, collider2
 double precision, intent(out) :: dratesq
 integer :: i,j,id1,p1,d1,id2,p2,d2
 double precision :: res !, delta
 res=0.0D0
! delta=0.0D0
 do i=1,peak1%Nc
  do j=1,peak2%Nc
   id1 = peak1%channel_id(i)
   p1 = int((id1-modulo(id1,10))/dble(10))
   d1 = modulo(id1,10)
   id2 = peak2%channel_id(j)
   p2 = int((id2-modulo(id2,10))/dble(10))
   d2 = modulo(id2,10)
!   if(p1.eq.p2) res=res+delta_rate%dCS(p1)**2*peak1%channel_w(i)*peak2%channel_w(j)
!   if(d1.eq.d2) res=res+delta_rate%dBR(d1)**2*peak1%channel_w(i)*peak2%channel_w(j)
   if((p1.gt.0).and.(p2.gt.0)) then
    if(collider1.ne.'ILC'.and.collider2.ne.'ILC') then
     if(delta_rate%CScov_ok.and.delta_rate%usecov) then
      if(useSMweights) then
        res=res+delta_rate%CScov(p1,p2)*peak1%channel_w(i)*peak2%channel_w(j)
       else 
        res=res+delta_rate%CScov(p1,p2)*peak1%channel_w_model(i)*peak2%channel_w_model(j)
      endif 
     else
      if(p1.eq.p2) then
       if(useSMweights) then
        res=res+delta_rate%dCS(p1)**2*peak1%channel_w(i)*peak2%channel_w(j)
       else
        res=res+delta_rate%dCS(p1)**2*peak1%channel_w_model(i)*peak2%channel_w_model(j)
       endif
      endif 
     endif
    else if(collider1.eq.'ILC'.and.collider2.eq.'ILC') then
     if(p1.eq.p2) then
      if(useSMweights) then
       res=res+delta_rate%dCS_ILC(p1)**2*peak1%channel_w(i)*peak2%channel_w(j)
      else
       res=res+delta_rate%dCS_ILC(p1)**2*peak1%channel_w_model(i)*peak2%channel_w_model(j)
      endif
     endif     
    endif 
   endif
!    if(collider1.ne.'ILC'.and.collider2.ne.'ILC') then
!     delta=delta_rate%dCS(p1)
!    else
!---Set the theoretical uncertainty of production mode to zero if collider==ILC     
!     delta=0.0D0
!    endif 
!    if(useSMweights) then
!     res=res+delta**2*peak1%channel_w(i)*peak2%channel_w(j)
!    else
!     res=res+delta**2*peak1%channel_w_model(i)*peak2%channel_w_model(j)
!    endif 
!   endif 
   if((d1.gt.0).and.(d2.gt.0)) then
    if(delta_rate%BRcov_ok.and.delta_rate%usecov) then
     if(useSMweights) then
       res=res+delta_rate%BRcov(d1,d2)*peak1%channel_w(i)*peak2%channel_w(j)
      else 
       res=res+delta_rate%BRcov(d1,d2)*peak1%channel_w_model(i)*peak2%channel_w_model(j)
     endif 
    else
     if(d1.eq.d2) then
      if(useSMweights) then
       res=res+delta_rate%dBR(d1)**2*peak1%channel_w(i)*peak2%channel_w(j)
      else 
       res=res+delta_rate%dBR(d1)**2*peak1%channel_w_model(i)*peak2%channel_w_model(j)
      endif 
     endif
    endif     
   endif   
  enddo
 enddo
  
 dratesq=res
end subroutine get_rate_uncertainties_sq_peaks


!------------------------------------------------------------------------------------       
subroutine get_cov_mu(matrix)
!------------------------------------------------------------------------------------       
 implicit none
 double precision, dimension(:,:), intent(out) :: matrix
 if(allocated(cov)) then
  if(size(cov,dim=1).eq.size(matrix,dim=1).and.size(cov,dim=2).eq.size(matrix,dim=2)) then
!  allocate(matrix(size(cov,dim=1),size(cov,dim=2)))
   matrix = cov
 else
  write(*,*) "WARNING in subroutine get_cov_mu: different dimensions."
  endif  
 else
  write(*,*) "WARNING in subroutine get_cov_mu: cov not allocated."
 endif
end subroutine get_cov_mu
!------------------------------------------------------------------------------------    
subroutine get_cov_mh(Hindex,matrix)
!------------------------------------------------------------------------------------       
 implicit none
 double precision, dimension(:,:), intent(out) :: matrix
 integer, intent(in) :: Hindex

 if(allocated(cov_mhneut)) then
  if(size(cov_mhneut,dim=2).eq.size(matrix,dim=1)&
&    .and.size(cov_mhneut,dim=3).eq.size(matrix,dim=2)) then
   if(Hindex.ge.lbound(cov_mhneut,dim=1).and.Hindex.le.ubound(cov_mhneut,dim=1)) then
    matrix = cov_mhneut(Hindex,:,:)
   else
    write(*,*) "WARNING in subroutine get_cov_mh: Hindex not in range."
   endif  
  else
   write(*,*) "WARNING in subroutine get_cov_mh: different dimensions."
  endif
 else
  write(*,*) "WARNING in subroutine get_cov_mh: cov not allocated."
 endif
 
end subroutine get_cov_mh
!------------------------------------------------------------------------------------       
subroutine print_cov_mu_to_file
 use usefulbits, only : file_id_common3
 implicit none
 character(LEN=20):: formatstring

 open(file_id_common3, file='cov_mu.txt', form='formatted')
 write(formatstring,'(A1,I2,A7)') '(',size(cov,dim=2),'F20.10)'
  do i=lbound(cov,dim=1),ubound(cov,dim=1)
   write(file_id_common3, formatstring) cov(i,:)
  enddo	
 close(file_id_common3)
end subroutine print_cov_mu_to_file
!------------------------------------------------------------------------------------    
subroutine print_corr_mu_to_file
 use usefulbits, only : file_id_common3
 implicit none
 character(LEN=20):: formatstring

 open(file_id_common3, file='corr_mu.txt', form='formatted')
 write(formatstring,'(A1,I2,A7)') '(',size(cov,dim=2),'F20.10)'
  do i=lbound(cov,dim=1),ubound(cov,dim=1)
! TS fix 14/08/2013: Took square root of this!  
   write(file_id_common3, formatstring) (cov(i,j)/sqrt(cov(i,i)*cov(j,j)), &
&   j=lbound(cov,dim=2),ubound(cov,dim=2))   
  enddo	
 close(file_id_common3)
 
end subroutine print_corr_mu_to_file
!------------------------------------------------------------------------------------    
subroutine print_inverse_cov_mu_to_file
 use usefulbits, only : file_id_common3
 implicit none
 character(LEN=20):: formatstring

 double precision, allocatable :: invcov(:,:)
 
 allocate(invcov(size(cov,dim=1),size(cov,dim=2)))
 call invmatrix(cov,invcov)  

 open(file_id_common3, file='inverse_cov_mu.txt', form='formatted')
 write(formatstring,'(A1,I2,A7)') '(',size(invcov,dim=2),'F20.10)'
  do i=lbound(invcov,dim=1),ubound(invcov,dim=1)
   write(file_id_common3, formatstring) invcov(i,:)
  enddo	
 close(file_id_common3)

 open(file_id_common3, file='mu_vector.txt', form='formatted')
 
  do i=lbound(mu_vector,dim=1),ubound(mu_vector,dim=1)
   write(file_id_common3, '(1F20.10)') mu_vector(i)
  enddo	
 close(file_id_common3)

 
end subroutine print_inverse_cov_mu_to_file
!------------------------------------------------------------------------------------    
subroutine print_inverse_cov_mh_to_file(Hindex)
 use usefulbits, only : file_id_common3
 implicit none
 character(LEN=20):: formatstring, filename
 integer, intent(in) :: Hindex

 double precision, allocatable :: invcov(:,:)
 
 allocate(invcov(size(cov_mhneut,dim=2),size(cov_mhneut,dim=3)))
 call invmatrix(cov_mhneut(Hindex,:,:),invcov)  

 write(filename,'(A14,I1,A4)') 'inverse_cov_mh',Hindex,'.txt'
 open(file_id_common3, file=trim(adjustl(filename)), form='formatted')
 write(formatstring,'(A1,I2,A7)') '(',size(invcov,dim=2),'F20.10)'
  do i=lbound(invcov,dim=1),ubound(invcov,dim=1)
   write(file_id_common3, formatstring) invcov(i,:)
  enddo	
 close(file_id_common3)
 
end subroutine print_inverse_cov_mh_to_file
!------------------------------------------------------------------------------------    
subroutine print_cov_mh_to_file(Hindex)
 use usefulbits, only : file_id_common3
 implicit none
 character(LEN=20):: formatstring, filename
 integer, intent(in) :: Hindex

 write(filename,'(A6,I1,A4)') 'cov_mh',Hindex,'.txt'
 open(file_id_common3, file=trim(adjustl(filename)), form='formatted')
 write(formatstring,'(A1,I2,A7)') '(',size(cov_mhneut,dim=3),'F20.10)'
  do i=lbound(cov_mhneut,dim=2),ubound(cov_mhneut,dim=2)
   write(file_id_common3, formatstring) cov_mhneut(Hindex,i,:)
  enddo	
 close(file_id_common3)
 
end subroutine print_cov_mh_to_file
!------------------------------------------------------------------------------------    
subroutine print_peakinformation
!------------------------------------------------------------------------------------    
! use usefulbits_HS, only : peaklist
 implicit none
 
 do i=lbound(analyses,dim=1),ubound(analyses,dim=1)
  do j=lbound(analyses(i)%peaks,dim=1),ubound(analyses(i)%peaks,dim=1)
   write(*,*)
   write(*,*)'#*************************************************************************#'   
   write(*,'(A,I3,A,I3,A)') ' #                        Analysis ',i,					&
&							'   Peak ',j,'                          #'
   write(*,*)'#*************************************************************************#'      
   write(*,'(A25,1I10)') 'ID =', analyses(i)%id
   write(*,'(A25,4X,A3)') 'Collaboration =', analyses(i)%table%collaboration
   write(*,'(A25,1F10.2)') 'Energy =', analyses(i)%table%energy   
   write(*,'(A25,4X,A45)') 'Reference =', analyses(i)%table%label
   write(*,'(A25,4X,A100)') 'Description =', analyses(i)%table%desc      
   write(*,'(A25,1F10.3)') 'mass resolution =',analyses(i)%table%deltam   
   write(*,'(A25,1F10.3)') 'peak mass =',analyses(i)%peaks(j)%mpeak
   write(*,'(A25,1F10.4)') 'peak mu =',analyses(i)%peaks(j)%mu
   write(*,'(A25,2F10.4)') 'cyan band(low,high) =',analyses(i)%peaks(j)%dmulow, 		&
&						   analyses(i)%peaks(j)%dmuup
!   write(*,'(A25,2F10.4)') 'dmu0 (low,high) =',analyses(i)%peaks(j)%dmulow0,			&
!&						   analyses(i)%peaks(j)%dmuup0
   write(*,'(A25,4X,$)')   'Higgs combination ='
   do k=lbound(analyses(i)%peaks(j)%Higgs_comb,dim=1),									&
&		ubound(analyses(i)%peaks(j)%Higgs_comb,dim=1)
    write(*,'(1I3,$)')		analyses(i)%peaks(j)%Higgs_comb(k)
   enddo
   write(*,*) 
   write(*,'(A25,7X,1I3)') 'Dominant Higgs =',analyses(i)%peaks(j)%domH
   write(*,'(A25,1F10.4)'), 'Total pred. mu =',analyses(i)%peaks(j)%total_mu
   write(*,'(A25,1F15.9)'), 'Chisq for mu =',analyses(i)%peaks(j)%chisq_mu
   write(*,'(A25,1F15.9)'), 'Chisq for mh =',analyses(i)%peaks(j)%chisq_mh
   write(*,'(A25,1F15.9)'), 'Chisq (total) =',analyses(i)%peaks(j)%chisq_tot
   write(*,'(A25,1F15.9)'), 'Chisq (max) =',analyses(i)%peaks(j)%chisq_max
   write(*,*)'#------------------------ Channel information ----------------------------#'
   write(*,*) '  ID       prod.	decay 	       mu  	     weight  	 syst.err. '
   write(*,*)'#-------------------------------------------------------------------------#'   
   do k=1, analyses(i)%peaks(j)%Nc
    write(*,'(1I5,5X,2A,3F15.6)') analyses(i)%peaks(j)%channel_id(k),&
    &	analyses(i)%table%channel_description(k,:),&
    &	analyses(i)%peaks(j)%channel_mu(k),analyses(i)%peaks(j)%channel_w_model(k),&
    &	analyses(i)%peaks(j)%channel_syst(k)
   enddo 
   write(*,*)'#-------------------------------------------------------------------------#'   
  enddo
 enddo 
end subroutine print_peakinformation  
!------------------------------------------------------------------------------------    
subroutine print_peakinformation_essentials
!------------------------------------------------------------------------------------    
!x use usefulbits_HS, only : peaklist
 implicit none
 
 do i=lbound(analyses,dim=1),ubound(analyses,dim=1)
  do j=lbound(analyses(i)%peaks,dim=1),ubound(analyses(i)%peaks,dim=1)
   write(*,*)
   write(*,*)'#*************************************************************************#'   
   write(*,'(A,I3,A,I3,A)') ' #                        Analysis ',i,					&
   &	 					'   Peak ',j,'                          #'
   write(*,*)'#*************************************************************************#'      
   write(*,'(A25,1I10)') 'ID =', analyses(i)%id
   write(*,'(A25,7X,A3)') 'Collaboration =', analyses(i)%table%collaboration   
   write(*,'(A25,4X,F6.2)') 'cms energy =', analyses(i)%table%energy   
   write(*,'(A25,4X,A45)') 'Reference =', analyses(i)%table%label
   write(*,'(A25,4X,A100)') 'Description =', analyses(i)%table%desc
   write(*,'(A25,1F10.2)') 'mass resolution =',analyses(i)%table%deltam         
   write(*,'(A25,1F10.2)') 'peak mass =',analyses(i)%peaks(j)%mpeak
   write(*,'(A25,1F10.4)') 'peak mu =',analyses(i)%peaks(j)%mu
   write(*,'(A25,2F10.4)') 'cyan band(low,high) =',analyses(i)%peaks(j)%dmulow, 		&
&							analyses(i)%peaks(j)%dmuup
   write(*,*)'#------------------------ Channel information ----------------------------#'
   write(*,*)'  ID       prod.	decay 	    efficiency'
   write(*,*)'#-------------------------------------------------------------------------#'   
   do k=1, analyses(i)%peaks(j)%Nc
    write(*,'(1I5,5X,2A,1F15.6)') analyses(i)%peaks(j)%channel_id(k),					&
&		analyses(i)%table%channel_description(k,:),analyses(i)%table%channel_eff(k)
   enddo 
   write(*,*)'#-------------------------------------------------------------------------#'   
  enddo
 enddo 
end subroutine print_peakinformation_essentials
!------------------------------------------------------------------------------------    
subroutine print_peaks_to_file
!------------------------------------------------------------------------------------    
 use usefulbits, only : file_id_common3
 use usefulbits_hs, only : StrCompress
 implicit none
 character(LEN=100) :: formatspec
 integer :: kk
 kk=0
 formatspec='(I3,7X,I10,1X,F6.2,1X,4F8.4,1X,A3,1X,F6.2,1X,F6.2,1X,I3,1X,A,5X,A)'
 open(file_id_common3,file="peak_information.txt")
 write(file_id_common3,*) "#HiggsSignals-"//trim(adjustl(HSvers))//						&
&						  " with experimental dataset '"//trim(adjustl(Exptdir))//"'" 
 write(file_id_common3,*) "#Number Analysis-ID mh_obs   mu_obs dmu_low dmu_high ",		&
 &				"dmh_exp  collaboration   energy	luminosity   description   reference"
 write(file_id_common3,*) "#"
 do i=lbound(analyses,dim=1),ubound(analyses,dim=1)
  do j=lbound(analyses(i)%peaks,dim=1),ubound(analyses(i)%peaks,dim=1)
  kk=kk+1
  write(file_id_common3,formatspec)kk,analyses(i)%id,analyses(i)%peaks(j)%mpeak,		&
  & analyses(i)%peaks(j)%mu, analyses(i)%peaks(j)%dmulow,analyses(i)%peaks(j)%dmuup,	&
  & analyses(i)%table%deltam,analyses(i)%table%collaboration, analyses(i)%table%energy,	&
  & analyses(i)%table%lumi, analyses(i)%table%mhchisq, 									&
  & trim(strcompress(analyses(i)%table%desc)), analyses(i)%table%label
  enddo
 enddo 
 close(file_id_common3)
end subroutine print_peaks_to_file
!------------------------------------------------------------------------------------    
subroutine print_peaks_to_LaTeX
!------------------------------------------------------------------------------------    
 use usefulbits, only : file_id_common3
 use usefulbits_hs, only : StrCompress
 implicit none
 character(LEN=100) :: formatspec
 integer :: kk, N, ii, id, p
 double precision :: weights(5) = (/ 0.0D0, 0.0D0,0.0D0,0.0D0,0.0D0 /)
 kk=0
 open(file_id_common3,file="peak_information.tex")
 write(file_id_common3,*) "\begin{tabular}{lcrrrrr}"
 write(file_id_common3,*) "\hline"
 write(file_id_common3,*) "Analysis & Signal strength & \multicolumn{5}{c}{Signal contamination [in \%]} \\"
 write(file_id_common3,*) "& & ggH & VBF & WH &  ZH & $t\bar{t}H$ \\"
 write(file_id_common3,*) "\hline"
 do i=lbound(analyses,dim=1),ubound(analyses,dim=1)
  do j=lbound(analyses(i)%peaks,dim=1),ubound(analyses(i)%peaks,dim=1)
  kk=kk+1
  N = analyses(i)%peaks(j)%Nc
  weights = (/ 0.0D0, 0.0D0,0.0D0,0.0D0,0.0D0 /)
  do ii=1,N
   id = analyses(i)%peaks(j)%channel_id(ii)
   p = int((id-modulo(id,10))/dble(10))
   weights(p) = analyses(i)%peaks(j)%channel_w_model(ii)
  enddo
  write(formatspec,"(A,I1,A,I1,A)") '(A3,1X,A,A,A,A,A,F6.2,A,F6.2,A,F6.2,A,F6.1,A,F6.1,A,F6.1,A,F6.1,A,F6.1,A,',N,'I3)'
  write(file_id_common3,formatspec) analyses(i)%table%collaboration, &
  & "$",trim(strcompress(analyses(i)%table%desc)), "$~\cite{", &
  & trim(strcompress(analyses(i)%table%label)),"} & $", &
  & analyses(i)%peaks(j)%mu, "\substack{+",analyses(i)%peaks(j)%dmuup,"\\ -",&
  & abs(analyses(i)%peaks(j)%dmulow),"}$ & $",100*weights(1),"$ & $",100*weights(2),&
  & "$ & $",100*weights(3),"$ & $",100*weights(4),"$ & $",100*weights(5),"$\\ %", &
  & analyses(i)%peaks(j)%channel_id
  enddo
 enddo 
 write(file_id_common3,*) "\hline"
 write(file_id_common3,*) "\end{tabular}"

end subroutine print_peaks_to_LaTeX
!------------------------------------------------------------------------------------    
subroutine print_peaks_signal_rates_to_file
!------------------------------------------------------------------------------------    
 use usefulbits, only : file_id_common3 
 use usefulbits_hs, only : HSres
 implicit none
 character(LEN=100) :: formatspec,formatspec2
 integer :: kk
 double precision :: mh_pull, mu_pull, dmu
 kk=0
 formatspec='(I3,7X,I10,1X,4F8.2,1X,6F10.4)'
 formatspec2='(I3,7X,I10,1X,2F8.2,1X,A7,1X,A7,1X,A7,1X,5F10.4)'

 open(file_id_common3,file="peak_massesandrates.txt")
 write(file_id_common3,*) "#HiggsSignals-"//trim(adjustl(HSvers))//						&
&						  " with experimental dataset '"//trim(adjustl(Exptdir))//"'" 
 write(file_id_common3,*) "#pull = (predicted - observed)/(gaussian uncertainty)"
 write(file_id_common3,*) "#Number Analysis-ID mh_obs dmh_exp mh_pred dmh_theo mh_pull",&
&						  " mu_obs dmu_low dmu_high mu_pred mu_pull"
 write(file_id_common3,*) "#"
 do i=lbound(analyses,dim=1),ubound(analyses,dim=1)
  do j=lbound(analyses(i)%peaks,dim=1),ubound(analyses(i)%peaks,dim=1)
  kk=kk+1
  if(analyses(i)%peaks(j)%domH.ne.0) then
   mh_pull = (analyses(i)%peaks(j)%Higgses(analyses(i)%peaks(j)%domH)%m - 				&
&			 analyses(i)%peaks(j)%mpeak)/												&
&			 sqrt(analyses(i)%table%deltam**2+											&
&			 analyses(i)%peaks(j)%Higgses(analyses(i)%peaks(j)%domH)%dm**2)
  endif
  call get_dmu_peak(dmu,analyses(i)%peaks(j))
  mu_pull = (analyses(i)%peaks(j)%total_mu - analyses(i)%peaks(j)%mu)/dmu
  
  if(analyses(i)%peaks(j)%domH.ne.0) then
  write(file_id_common3,formatspec)kk,analyses(i)%id,analyses(i)%peaks(j)%mpeak,		&
&  analyses(i)%table%deltam, analyses(i)%peaks(j)%Higgses(analyses(i)%peaks(j)%domH)%m,	&
&  analyses(i)%peaks(j)%Higgses(analyses(i)%peaks(j)%domH)%dm, mh_pull,					&
&  analyses(i)%peaks(j)%mu, analyses(i)%peaks(j)%dmulow,analyses(i)%peaks(j)%dmuup,		&
&  analyses(i)%peaks(j)%total_mu, mu_pull
  else	
    write(file_id_common3,formatspec2)kk,analyses(i)%id,analyses(i)%peaks(j)%mpeak,		&
&  analyses(i)%table%deltam, 'NAN','NAN','NAN',											&
&  analyses(i)%peaks(j)%mu, analyses(i)%peaks(j)%dmulow,analyses(i)%peaks(j)%dmuup,		&
&  analyses(i)%peaks(j)%total_mu, mu_pull
  endif	
  enddo
 enddo 
 close(file_id_common3)
end subroutine print_peaks_signal_rates_to_file
!--------------------------------------------------------------------
subroutine get_peakchi2(obsID, csqmu, csqmh, csqmax, csqtot)
!--------------------------------------------------------------------
 implicit none
 integer, intent(in) :: obsID
 double precision, intent(out) :: csqmu, csqmh, csqmax, csqtot
 integer :: i, j
 
 do i=lbound(analyses,dim=1),ubound(analyses,dim=1)
  do j=lbound(analyses(i)%peaks,dim=1),ubound(analyses(i)%peaks,dim=1)
   if(obsID.eq.analyses(i)%peaks(j)%id) then
    csqmu = analyses(i)%peaks(j)%chisq_mu
    csqmh = analyses(i)%peaks(j)%chisq_mh
    csqmax = analyses(i)%peaks(j)%chisq_max
    csqtot = analyses(i)%peaks(j)%chisq_tot
   endif
  enddo
 enddo

end subroutine get_peakchi2
!------------------------------------------------------------------------------------     
subroutine add_peaks_to_HSresults(r)
!------------------------------------------------------------------------------------    
! use usefulbits_hs, only : HSresults, peaklist
 implicit none
 
 type(HSresults), intent(out) :: r
 
 integer :: i,j, iii
 if(allocated(r%obsID)) deallocate(r%obsID)
 if(allocated(r%mupred)) deallocate(r%mupred)
 if(allocated(r%domH)) deallocate(r%domH)
 if(allocated(r%nH)) deallocate(r%nH)

!-The number of peaks should maybe be saved in a global usefulbits_HS variable in order
! to avoid multiple loops over the peaklist.
 iii=0
 do i=lbound(analyses,dim=1),ubound(analyses,dim=1)
  do j=lbound(analyses(i)%peaks,dim=1),ubound(analyses(i)%peaks,dim=1)
   iii=iii+1
  enddo
 enddo
 
 allocate(r%mupred(iii), r%domH(iii), r%nH(iii), r%obsID(iii)) 
  
 iii=0
 do i=lbound(analyses,dim=1),ubound(analyses,dim=1)
  do j=lbound(analyses(i)%peaks,dim=1),ubound(analyses(i)%peaks,dim=1)
   iii=iii+1
   r%mupred(iii)=analyses(i)%peaks(j)%total_mu
   r%domH(iii)=analyses(i)%peaks(j)%domH
   r%nH(iii)=analyses(i)%peaks(j)%NHiggs_comb
   r%obsID(iii)=analyses(i)%peaks(j)%id
  enddo
 enddo
 
end subroutine add_peaks_to_HSresults
!------------------------------------------------------------------------------------    
subroutine get_dmu0sq_peak(dmu0sq,peak,domax)
!------------------------------------------------------------------------------------    
 use usefulbits_HS, only : mupeak, symmetricerrors, absolute_errors
 implicit none
 
 integer, intent(in) :: domax   ! if 1, then use predicted mu == 0 
 double precision, intent(out) :: dmu0sq
 double precision :: pred_mu, mu
 
 type(mupeak), intent(in) :: peak

! TESTING:
! if(absolute_errors) then
!  mu=peak%mu_original
! else
 mu=peak%mu
! endif

 if(domax.eq.1) then
  pred_mu = 0.0D0
 else
  pred_mu = peak%total_mu 
 endif
 
 if(.not.symmetricerrors) then
  if(pred_mu.le.mu) then
   dmu0sq = peak%dmulow0sq
  else if(pred_mu.gt.mu) then
   dmu0sq = peak%dmuup0sq
  endif 
 else
  if(peak%dmulow0sq.lt.0.0d0.or.peak%dmuup0sq.lt.0.0d0) then
   write(*,*) "WARNING: squared intrinsic mu uncertainty is negative!"
   dmu0sq= (sqrt(abs(peak%dmulow0sq))+sqrt(abs(peak%dmuup0sq)))**2/4.
  else 
  dmu0sq = (sqrt(peak%dmulow0sq)+sqrt(peak%dmuup0sq))**2/4.
  endif
 endif
 
!! write(*,*) 'DEBUG: ',peak%dmulow0sq, peak%dmuup0sq, (sqrt(peak%dmulow0sq)+sqrt(peak%dmuup0sq))**2/4.
 
end subroutine get_dmu0sq_peak
!------------------------------------------------------------------------------------    
subroutine get_dmu_peak(dmu,peak)
!------------------------------------------------------------------------------------    
 use usefulbits_HS, only : mupeak
 implicit none
 
 double precision, intent(out) :: dmu
 double precision :: pred_mu
 
 type(mupeak), intent(in) :: peak
 
 pred_mu = peak%total_mu 
 
 if(pred_mu.le.peak%mu) then
  dmu = peak%dmulow
 else if(pred_mu.gt.peak%mu) then
  dmu = peak%dmuup
 endif 

end subroutine get_dmu_peak
!------------------------------------------------------------------------------------    
subroutine calc_pc_chisq(pc_chisq,peak,mhchisqflag,Higgses,Higgs_dom,indices_in,iterstep)
! This subroutine calculates the chisq value for a given Higgs combination assigned
! to one peak.
! Parameters:
!
!  pc_chisq:       Return value (double precision)
!  peak:           Peak observable the Higgs bosons are assigned to (type mupeak)
!  Higgses:		   neutral Higgs bosons considered for the assignment (type neutHiggs(:))
!  Hdom:		   dominantly contribution Higgs boson (int)	
!  indices_in:     Higgs boson combination the chi squared value is evaluated for (int(:))
!  iterstep:	   Iteration-step of Higgs-to-peak assignment
!------------------------------------------------------------------------------------    
 implicit none
 type(mupeak), intent(in) :: peak
 type(neutHiggs), dimension(:), intent(in) :: Higgses(:)
! integer, allocatable, intent(out) :: indices_best(:)
! integer, intent(in), optional :: Nindices_in
 integer, dimension(:), intent(in) :: indices_in(:)
 integer, intent(in) :: mhchisqflag, iterstep
 integer, intent(out) :: Higgs_dom
 double precision, intent(out) :: pc_chisq
 
!--Internal parameters:
 double precision :: mutot,dmu,csq0,csq_mh_tot,csq_tot,csq_tot_tmp,csq_mh_tmp,mumax
 double precision :: deltam
 integer :: nH
 
 call check_pdf(pdf)
 
 nH=size(Higgses)
!--Calculate the maximal chisq value for this peak
 csq0 = csq_mu(0.0D0,peak%mu,peak%dmuup,peak%dmulow)

 mumax = -1.0D6
 mutot = 0.0D0
 csq_mh_tot=0.0D0
 Higgs_dom=0
  
 do i=lbound(indices_in,dim=1),ubound(indices_in,dim=1)
  if(indices_in(i).ne.0) then  
   mutot = mutot + Higgses(indices_in(i))%mu
!--Determining the dominantly contributing Higgs boson. 
   if(Higgses(indices_in(i))%mu.gt.mumax) then
    mumax = Higgses(indices_in(i))%mu
    Higgs_dom = indices_in(i)
   endif
   if(mhchisqflag.eq.1) then
    deltam = peak%dm
    csq_mh_tmp = csq_mh(Higgses(indices_in(i))%m, peak%mpeak,							&
& 					     Higgses(indices_in(i))%dm,deltam)
   else
    csq_mh_tmp = 0.0D0
   endif   
   csq_mh_tot = csq_mh_tot + csq_mh_tmp
   
  endif	
 enddo

 if(allocated(cov_mhneut_max).and.iterstep.eq.1) then
 !-In the first iterated step, use a Higgs mass covariance matrix, where all Higgs bosons
 !-are assumed to be assigned.
  csq_mh_tot = csq_mh_with_max_corr(indices_in, peak%internalnumber)
 elseif(allocated(cov_mhneut).and.iterstep.gt.1) then
!-In the second (and succeeding) iterated steps, use a the Higgs mass covariance matrix
!-based on the the previous Higgs-to-peaks assignment.
  csq_mh_tot = csq_mh_with_corr(indices_in, peak%internalnumber)  
 endif
 if(allocated(cov)) then
  csq_tot = csq_mh_tot + csq_mu_with_corr(mutot, peak%internalnumber)
 else
  csq_tot = csq_mh_tot + csq_mu(mutot, peak%mu, peak%dmuup, peak%dmulow) 
 endif

 pc_chisq = csq_tot  
  
end subroutine calc_pc_chisq
!---------------------------------------------------------------------------
function csq_mh_with_max_corr(indices_in, peaknumber)
!---------------------------------------------------------------------------
 integer, dimension(:), intent(in) :: indices_in(:)
 integer, intent(in) :: peaknumber
 double precision :: csq_mh_with_max_corr 
 double precision, allocatable :: csq_mh_per_Higgs(:)
 integer :: i, ii, iii, k, nH, N, Hindex
 double precision, allocatable :: v(:,:), invcov(:,:), v2(:)

  nH = size(cov_mhneut,dim=1)
  N = size(cov_mhneut,dim=2)
  allocate(v(nH,N), v2(N), csq_mh_per_Higgs(nH))
  
  do k=1,nH
   iii=0
   do i=1, size(analyses)
    do ii=lbound(analyses(i)%peaks,dim=1),ubound(analyses(i)%peaks,dim=1)
     iii=iii+1
      v(k,iii) = analyses(i)%peaks(ii)%Higgses(k)%m - analyses(i)%peaks(ii)%mpeak
      if(iii.eq.peaknumber) v(k,iii) = 0.0D0 ! Will be filled later...
    enddo
   enddo
  enddo 

  do k=1,nH
   iii=0
   do i=1, size(analyses)
    do ii=lbound(analyses(i)%peaks,dim=1),ubound(analyses(i)%peaks,dim=1)
     iii=iii+1
     if(iii.eq.peaknumber) then
      Hindex=indices_in(k)
      if(Hindex.ne.0) then     
       v(Hindex,iii) = analyses(i)%peaks(ii)%Higgses(k)%m - analyses(i)%peaks(ii)%mpeak
      endif  
     endif 
    enddo
   enddo
  enddo

   do k=1,nH 		!-n.b.: this loops now over Hindex
   call invmatrix(cov_mhneut_max(k,:,:),invcov) 
   call matmult(invcov,v(k,:),v2,N,1)
   csq_mh_per_Higgs(k) = v(k,peaknumber)*v2(peaknumber)
  enddo  
  
csq_mh_with_max_corr = sum(csq_mh_per_Higgs)

deallocate(v,v2,csq_mh_per_Higgs)

end function csq_mh_with_max_corr
!---------------------------------------------------------------------------
function csq_mh_with_corr(indices_in, peaknumber)
!---------------------------------------------------------------------------
 integer, dimension(:), intent(in) :: indices_in(:)
 integer, intent(in) :: peaknumber
 double precision :: csq_mh_with_corr 
 double precision, allocatable :: csq_mh_per_Higgs(:)
 integer :: i, ii, iii, k, nH, N, Hindex
 double precision, allocatable :: v(:,:), invcov(:,:), v2(:)

 nH = size(cov_mhneut,dim=1)
 N = size(cov_mhneut,dim=2)
 allocate(v(nH,N), v2(N), csq_mh_per_Higgs(nH))
    
!-First, fill the vectors with zeros:
 do k=1,nH
  do i=1,N
   v(k,i) = 0.0D0
  enddo
 enddo  

  do k=1,nH
   iii=0
   do i=1, size(analyses)
    do ii=lbound(analyses(i)%peaks,dim=1),ubound(analyses(i)%peaks,dim=1)
     iii=iii+1
     Hindex = analyses(i)%peaks(ii)%Higgs_comb(k)
     if(iii.eq.peaknumber) Hindex=indices_in(k)     
     if(Hindex.ne.0) then     
      v(Hindex,iii) = analyses(i)%peaks(ii)%Higgses(Hindex)%m -							&
&					  analyses(i)%peaks(ii)%mpeak
     endif  
    enddo
   enddo
  enddo 

  do k=1,nH 		!-n.b.: this loops now over Hindex
   call invmatrix(cov_mhneut(k,:,:),invcov) 
   call matmult(invcov,v(k,:),v2,N,1)
   csq_mh_per_Higgs(k) = v(k,peaknumber)*v2(peaknumber)
  enddo  
  
  csq_mh_with_corr = sum(csq_mh_per_Higgs)

  deallocate(v,v2,csq_mh_per_Higgs)

 end function csq_mh_with_corr
!------------------------------------------------------------------------------------    
function csq_mu_with_corr(mu, peaknumber)
!------------------------------------------------------------------------------------    
! use usefulbits_hs, only : peaklist, cov
 use numerics, only : invmatrix, matmult

 integer, intent(in) :: peaknumber
 double precision, intent(in) :: mu
 integer :: i, ii, iii, N
 double precision, allocatable :: v(:), vmat(:,:), invcov(:,:), v2(:)
 double precision :: csq_mu_with_corr
 
 if(allocated(cov)) deallocate(cov)
 call create_covariance_matrix_mu(0) 
 
 N = size(cov,dim=1)
 allocate(v(N), vmat(N,1),invcov(N,N), v2(N))
 
 !-First construct the vector (mupred - muobs)_iii
 iii=0
 do i=1, size(analyses)
  do ii=lbound(analyses(i)%peaks,dim=1),ubound(analyses(i)%peaks,dim=1)
   iii=iii+1
   v(iii) = analyses(i)%peaks(ii)%total_mu - analyses(i)%peaks(ii)%mu
   if(iii.eq.peaknumber) v(iii) = mu - analyses(i)%peaks(ii)%mu
   vmat(iii,1) = v(iii)
  enddo
 enddo
    
 call invmatrix(cov,invcov)   
 call matmult(invcov,vmat,v2,N,1)
 
 csq_mu_with_corr = v(peaknumber)*v2(peaknumber)

 deallocate(v,vmat,invcov,v2)

end function csq_mu_with_corr
!------------------------------------------------------------------------------------  
subroutine check_pdf(pdf)
!------------------------------------------------------------------------------------  
 implicit none
 integer, intent(inout) :: pdf
 
 if(.not.((pdf.eq.1).or.(pdf.eq.2).or.(pdf.eq.3))) then
  write(*,*) 'WARNING: pdf not properly specified. Will be set to pdf=2 (gaussian-shape)'
  pdf=2
 endif
  
end subroutine check_pdf
!------------------------------------------------------------------------------------    
function csq_mu(x, x0, dx_up, dx_low)
!- x: model predicted value
!- x0: observed value
!- dx_up: difference between observed and upper 1sigma band (always positive)
!- dx_low: difference between observed and lower 1sigma band (always positive)
!------------------------------------------------------------------------------------    
 use usefulbits_hs, only : symmetricerrors
 implicit none

 double precision, intent(in) :: x, x0, dx_up, dx_low
 double precision :: csq_mu, dx
 
 if(.not.symmetricerrors) then
  if(x.ge.x0) then
   dx = dx_up
  else
   dx = dx_low
  endif
 else
  dx = (abs(dx_up)+abs(dx_low))/2.
 endif
  
 csq_mu=chisq(x, x0, dx)
 
end function csq_mu 
!------------------------------------------------------------------------------------      
function csq_mh(x,x0,dx1,dx2)
!------------------------------------------------------------------------------------    
 implicit none
 double precision, intent(in) :: x, x0, dx1, dx2
 double precision :: csq_mh, dx
 
 if(pdf.eq.1) then			! box pdf
  dx = dx1 + dx2			! Add exp. and th. uncertainties linearly
  if(x.ge.(x0-dx).and.x.le.(x0+dx)) then
   csq_mh = 0.0D0
  else
   csq_mh = 100000.0D0	!Large number
  endif
 else if(pdf.eq.2) then		! gaussian pdf
  dx = sqrt(dx1**2+dx2**2)	! Add uncertainties in quadrature
  csq_mh = chisq(x,x0,dx)
!---------------------------- Have to verify boxgaussian...
 else if(pdf.eq.3) then		! box+gaussian pdf
  if(x.ge.(x0-dx1).and.x.le.(x0+dx1)) then
   csq_mh = 0.0D0
  else if(x.ge.(x0+dx1)) then
   csq_mh = chisq(x,x0+dx1,dx2)
  else if(x.le.(x0-dx1)) then
   csq_mh = chisq(x,x0-dx1,dx2)
  endif
 endif 
end function csq_mh 
!------------------------------------------------------------------------------------    
function chisq(x,x0,dx)
!------------------------------------------------------------------------------------    
 implicit none 
 double precision, intent(in) :: x,x0,dx
 double precision :: chisq
 
 if(dx.ne.0) then 
  chisq=(x-x0)**2/dx**2
 else
  write(*,*) 'WARNING, dx = 0'
  chisq=100000. 
 endif
 
end function chisq
!------------------------------------------------------------------------------------      
subroutine get_ncomb(ncomb, indices_truncated, indices_best)
!------------------------------------------------------------------------------------     
 implicit none
 
 integer, dimension(:), intent(in) :: indices_best(:)
 integer, allocatable, intent(out) :: indices_truncated(:)
 integer, intent(out) :: ncomb
 integer :: i,j

 ncomb=0

 do i=lbound(indices_best,dim=1),ubound(indices_best,dim=1)
  if(indices_best(i).ne.0) ncomb=ncomb+1
 enddo
 
 allocate(indices_truncated(ncomb))
 j=1
 do i=lbound(indices_best,dim=1),ubound(indices_best,dim=1)
  if(indices_best(i).ne.0) then
   indices_truncated(j)=indices_best(i)
   j=j+1
  endif 
 enddo 
end subroutine get_ncomb 
!------------------------------------------------------------------------------------  
subroutine get_weights_at_peak( peak, mutab )
! This subroutines fills the channels weights array of the peak object with the
! Standard Model weights obtained at the peaks position.
!------------------------------------------------------------------------------------  
 use usefulbits, only : div, small
! use theory_XS_SM_functions
 use usefulbits_HS, only : mutable

 type(mupeak), intent(inout) :: peak
 type(mutable), intent(in) :: mutab

 integer :: i, id, p, d
 double precision :: SMrate, mass
 double precision :: SMCS_lhc7_gg_H,SMCS_lhc7_bb_H,SMCS_lhc7_vbf_H,SMCS_lhc7_HW,	&
 &	SMCS_lhc7_HZ, SMCS_lhc7_ttH, SMCS_lhc8_gg_H,SMCS_lhc8_bb_H,SMCS_lhc8_vbf_H,		&
 &	SMCS_lhc8_HW, SMCS_lhc8_HZ, SMCS_lhc8_ttH, SMCS_tev_gg_H,SMCS_tev_bb_H,			&
 &  SMCS_tev_vbf_H, SMCS_tev_HW, SMCS_tev_HZ, SMCS_tev_ttH,SMBR_Hgamgam,SMBR_HWW,	&
 &  SMBR_HZZ, SMBR_Htautau, SMBR_Hbb,SMBR_HZgam,SMBR_Hcc, SMBR_Hmumu,SMBR_Hgg
 ! Check experiment and energy flag to choose the relevant dataset

 mass = peak%mpeak

 do i=1,mutab%Nc
  id = mutab%channel_id(i)
  p = int((id-modulo(id,10))/dble(10))
  d = modulo(id,10)
  
!--Do the production rate for the relevant experiment and cms-energy 
  if(mutab%collider.eq.'LHC') then
   if(abs(mutab%energy-7.0D0).le.small) then
    if(p.eq.1) then 
     SMrate=SMCS_lhc7_gg_H(mass)+SMCS_lhc7_bb_H(mass)
    else if(p.eq.2) then
     SMrate=SMCS_lhc7_vbf_H(mass)
    else if(p.eq.3) then
     SMrate=SMCS_lhc7_HW(mass)
    else if(p.eq.4) then
     SMrate=SMCS_lhc7_HZ(mass)
    else if(p.eq.5) then
     SMrate=SMCS_lhc7_ttH(mass)
    endif 
   else if(abs(mutab%energy-8.0D0).le.small) then
    if(p.eq.1) then 
     SMrate=SMCS_lhc8_gg_H(mass)+SMCS_lhc8_bb_H(mass)
    else if(p.eq.2) then
     SMrate=SMCS_lhc8_vbf_H(mass)
    else if(p.eq.3) then
     SMrate=SMCS_lhc8_HW(mass)
    else if(p.eq.4) then
     SMrate=SMCS_lhc8_HZ(mass)
    else if(p.eq.5) then
     SMrate=SMCS_lhc8_ttH(mass)
    endif 
   endif
  else if(mutab%collider.eq.'TEV') then
    if(p.eq.1) then 
     SMrate=SMCS_tev_gg_H(mass)+SMCS_tev_bb_H(mass)
    else if(p.eq.2) then
     SMrate=SMCS_tev_vbf_H(mass)
    else if(p.eq.3) then
     SMrate=SMCS_tev_HW(mass)
    else if(p.eq.4) then
     SMrate=SMCS_tev_HZ(mass)
    else if(p.eq.5) then
     SMrate=SMCS_tev_ttH(mass)
    endif 
  endif
!--Multiply now by the decay rate
  if(d.eq.1) then
   SMrate=SMrate*SMBR_Hgamgam(mass)
  else if(d.eq.2) then
   SMrate=SMrate*SMBR_HWW(mass)
  else if(d.eq.3) then
   SMrate=SMrate*SMBR_HZZ(mass)
  else if(d.eq.4) then
   SMrate=SMrate*SMBR_Htautau(mass)
  else if(d.eq.5) then
   SMrate=SMrate*SMBR_Hbb(mass)
  else if(d.eq.6) then
   SMrate=SMrate*SMBR_HZgam(mass)
  else if(d.eq.7) then
   SMrate=SMrate*SMBR_Hcc(mass)
  else if(d.eq.8) then
   SMrate=SMrate*SMBR_Hmumu(mass)
  else if(d.eq.9) then
   SMrate=SMrate*SMBR_Hgg(mass)   
  endif
  
  peak%channel_w(i)=peak%channel_eff(i)*SMrate
 enddo
 
 SMrate=sum(peak%channel_w(:))
 do i=1,mutab%Nc
  peak%channel_w(i)=div(peak%channel_w(i),SMrate,0.0D0,1.0D9)
 enddo
 
end subroutine get_weights_at_peak
!------------------------------------------------------------------------------------   
end module pc_chisq
!------------------------------------------------------------------------------------  