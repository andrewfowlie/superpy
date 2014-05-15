!--------------------------------------------------------------------
! This file is part of HiggsSignals (OS and TS 04/03/2013)
!--------------------------------------------------------------------
module datatables
!--------------------------------------------------------------------
 use usefulbits_HS, only : mutable, mupeak, Exptdir, obs, analyses !, withcorrexpsyst
!x  mutables, peaks
 use usefulbits, only : np
 implicit none

 integer :: ntable1,ntable2

 contains
!--------------------------------------------------------------------
 subroutine initializemutable_blank(table)
!--------------------------------------------------------------------
  type(mutable) :: table

   table%id         = -1
   table%nx         = -1
   table%particle_x = -1
   table%label      = ''
   table%desc       = ''
   table%expt       = ''
   table%lumi       = -1.0D0
   table%dlumi		= -1.0D0
   table%energy     = -1.0D0
   table%Nc			= -1
   table%xmax       = -1.0D0
   table%xmin       = -1.0D0
   table%sep        = -1.0D0
   table%deltam     = -1.0D0
   table%deltax     = -1.0D0
   table%mhchisq     = 0
 
 end subroutine initializemutable_blank
!--------------------------------------------------------------------
 subroutine initialize_observables
!--------------------------------------------------------------------
  use store_pathname_HS
  use usefulbits, only: np,Hneut,Hplus,file_id_common3
  implicit none
  character(LEN=150) :: filename

  character(LEN=100) :: datafile(500)
  integer n_datafiles

  !--------------------------------------input
!  type(mutable),allocatable :: tables(:)
  !-----------------------------------internal
  logical :: newtables
  integer :: x, n, xbeg, xend, i
  character(LEN=pathname_length+150) :: fullfilename    
  character(LEN=200) :: comment
  character(LEN=1) :: firstchar
  character(LEN=100) :: line
  integer :: col
  integer :: ios
  integer, allocatable :: skip(:)
  !-------------------------------------------  
!  open(file_id_common3, file=trim(adjustl(pathname_HS))//							&
!  &							 "Expt_tables/analyses.txt", form='formatted')
!  read(file_id_common3,'(A)',iostat=ios) datafile(1)
	
!  if(ios.ne.0)then
!   close(file_id_common3)
!!   call system('ls -1 -p '//trim(adjustl(pathname_HS))//							&
!!&	 'Expt_tables/'//trim(adjustl(Exptdir))//' > '//trim(adjustl(pathname_HS))//	&
!!&	 'Expt_tables/analyses.txt')
   call system('ls -1 -p '//trim(adjustl(pathname_HS))// &
&	 'Expt_tables/'//trim(adjustl(Exptdir))//' > HS_analyses.txt')

!!   open(file_id_common3, file=trim(adjustl(pathname_HS))//"Expt_tables/analyses.txt",&
!!&    form='formatted')
   open(file_id_common3, file="HS_analyses.txt",form='formatted')

!  else
!   rewind(file_id_common3)
!  endif

  print *, "Reading in the following datafiles from analysis-set "//				&
  			trim(adjustl(Exptdir))//":"
  n = 0
  n_datafiles = 0
  do
   n = n+1
   read(file_id_common3,'(A)', iostat=ios) datafile(n)
   if(ios.ne.0) exit
   write(*,'(I4,2X,A)') n, datafile(n)
  enddo
  n_datafiles = n - 1

  close(file_id_common3)
	
  allocate(obs(n_datafiles),skip(n_datafiles))            
  do i=lbound(obs,dim=1), ubound(obs,dim=1)
   call initializemutable_blank(obs(i)%table)
  enddo

  do n=1,n_datafiles
   skip(n)=11
!!   if(withcorrexpsyst) skip(n)=12
   open(file_id_common3, file=trim(adjustl(pathname_HS)) //'Expt_tables/'// &
&	 trim(adjustl(Exptdir))//'/' // datafile(n))
   do 
   	read(file_id_common3,'(A)') comment
	comment = trim(adjustl(comment))
	write(firstchar,'(A1)') comment   
	if(firstchar.ne.'#') then
	 exit
    else
	 skip(n)=skip(n)+1
    endif  
   enddo 
   backspace(file_id_common3)
   read(file_id_common3,*) obs(n)%id, obs(n)%table%id, obs(n)%obstype
   read(file_id_common3,'(A)') obs(n)%table%label
   read(file_id_common3,*) obs(n)%table%collider,obs(n)%table%collaboration, &
   & obs(n)%table%expt
   read(file_id_common3,'(A)') obs(n)%table%desc
   read(file_id_common3,*) obs(n)%table%energy, obs(n)%table%lumi, obs(n)%table%dlumi
!--TESTING correlated experimental systematics:
!!   if(withcorrexpsyst) read(file_id_common3,*) (obs(n)%table%correxpsyst(i),i=1,4)
!!   write(*,*) "Systematics: ",obs(n)%table%correxpsyst
!--END   
   read(file_id_common3,*) obs(n)%table%particle_x, obs(n)%table%mhchisq
!--CHECK FOR ASSIGNMENT GROUP AS SECOND COLUMN:
   read(file_id_common3,'(A)') line
   call read_in_mass_resolution_and_assignment_group(line, obs(n)%table%deltam,&
&   obs(n)%table%assignmentgroup)
!   write(*,*) "dm, group = ",obs(n)%table%deltam, obs(n)%table%assignmentgroup
!   read(file_id_common3,*) obs(n)%table%deltam
   read(file_id_common3,*) obs(n)%table%xmin, obs(n)%table%xmax, obs(n)%table%sep
   read(file_id_common3,*) obs(n)%table%Nc, obs(n)%table%eff_ref_mass
   allocate(obs(n)%table%channel_id(obs(n)%table%Nc))
   read(file_id_common3,*) (obs(n)%table%channel_id(i),i=1,obs(n)%table%Nc)
   allocate(obs(n)%table%channel_eff(obs(n)%table%Nc))
   allocate(obs(n)%table%channel_eff_ratios(obs(n)%table%Nc))   
   if(obs(n)%table%eff_ref_mass.ge.0D0) then
    read(file_id_common3,*) (obs(n)%table%channel_eff(i),i=1,obs(n)%table%Nc)
   else
    do i=1,obs(n)%table%Nc
     obs(n)%table%channel_eff(i)=1.0D0
    enddo 
    read(file_id_common3,*)
   endif
   obs(n)%table%channel_eff_ratios=1.0D0
   allocate(obs(n)%table%channel_description(obs(n)%table%Nc,2))			
   close(file_id_common3)
  enddo
				
  col = 3
  do x=1,n_datafiles	
						
   if(obs(x)%table%deltam.le.obs(x)%table%sep) then
    write(*,*) "Warning: Mass resolution for " ,trim(obs(x)%table%label),			&
&    "observable number ",x," with ID = ",obs(x)%id,								&
&    " very small - may lead to unreliable results"
	write(*,*) "Mass resolution = ",obs(x)%table%deltam
	write(*,*) "Separation = ",obs(x)%table%sep
   endif

   obs(x)%table%nx=nint((obs(x)%table%xmax-obs(x)%table%xmin)/obs(x)%table%sep)+1
   allocate(obs(x)%table%mu(obs(x)%table%nx,col))  			
   allocate(obs(x)%table%mass(obs(x)%table%nx))   
   allocate(obs(x)%table%channel_mu(obs(x)%table%Nc,np(obs(x)%table%particle_x)))
   allocate(obs(x)%table%channel_systSM(obs(x)%table%Nc,np(obs(x)%table%particle_x)))
   allocate(obs(x)%table%channel_syst(obs(x)%table%Nc,np(obs(x)%table%particle_x)))
   allocate(obs(x)%table%channel_w(obs(x)%table%Nc,np(obs(x)%table%particle_x)))
   allocate(obs(x)%table%channel_w_corrected_eff(obs(x)%table%Nc,np(obs(x)%table%particle_x)))
  enddo
  		
   do x=1, n_datafiles
    fullfilename=trim(adjustl(pathname_HS))//'Expt_tables/'//trim(adjustl(Exptdir))//'/'&
&                //trim(datafile(x))
    call read_mutable(obs(x)%table,skip(x),col,fullfilename) 
   enddo
  deallocate(skip)
  		
 end subroutine initialize_observables
!--------------------------------------------------------------------
 subroutine read_in_mass_resolution_and_assignment_group(line, dm, group)
!--------------------------------------------------------------------
 character(LEN=100), intent(in) :: line
 character(LEN=100), intent(out) :: group
 double precision, intent(out) :: dm
 
 integer :: i, indx, prev, beginning
 integer :: j, indxstr, prevstr, beginningstr

 prev      = 0 
 beginning = 1 
 prevstr      = 0
 beginningstr = 1

 group = ""
 
 do i=1,len(line)

  indx = index('0123456789.', line(i:i))
  if (indx.eq.0 .and. prev.gt.0) then
   read(line(beginning:i-1), *) dm
  else if (indx.gt.0 .and. prev.eq.0) then
   beginning = i 
  end if
  prev = indx
         
  indxstr = index('ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz_', line(i:i))
  if (indxstr.eq.0 .and. prevstr.gt.0) then
   read(line(beginningstr:i-1), *) group
  else if (indxstr.gt.0 .and. prevstr.eq.0) then
   beginningstr = i
  end if
  prevstr = indxstr
 end do
 
! group=trim(adjustl(group_tmp))
 
 end subroutine read_in_mass_resolution_and_assignment_group 
!--------------------------------------------------------------------
 subroutine setup_peak_observables
!--------------------------------------------------------------------
  implicit none
  integer :: i,j
  
  do i=lbound(obs,dim=1),ubound(obs,dim=1)
   if(obs(i)%obstype.eq.1) then
    obs(i)%peak%id = obs(i)%id
    obs(i)%peak%mpeak = obs(i)%table%mass(1)
    obs(i)%peak%dm = obs(i)%table%deltam
    obs(i)%peak%mu = obs(i)%table%mu(1,2)
    obs(i)%peak%mu_original = obs(i)%table%mu(1,2)    
    obs(i)%peak%dmuup = obs(i)%table%mu(1,3) - obs(i)%table%mu(1,2)
    obs(i)%peak%dmulow = obs(i)%table%mu(1,2) - obs(i)%table%mu(1,1)
    obs(i)%peak%dlumi = obs(i)%table%dlumi
	obs(i)%peak%scale_mu = 1.0D0
	obs(i)%peak%undo_assignment = 0
	obs(i)%peak%assignmentgroup = obs(i)%table%assignmentgroup
!	obs(i)%peak%scale_mh = 1.0D0	
    if(.not.allocated(obs(i)%peak%Higgs_comb)) then
     allocate(obs(i)%peak%Higgs_comb(np(obs(i)%table%particle_x)))
    endif 
	if(.not.allocated(obs(i)%peak%Higgses)) then
	 allocate(obs(i)%peak%Higgses(np(obs(i)%table%particle_x)))
	endif 
	obs(i)%peak%Nc = obs(i)%table%Nc
   ! have to allocate channel info and set them initially to zero!
	if(.not.allocated(obs(i)%peak%channel_id)) then
	 allocate(obs(i)%peak%channel_id(obs(i)%peak%Nc))
	endif 
    obs(i)%peak%channel_id(:) = obs(i)%table%channel_id(:)
	if(.not.allocated(obs(i)%peak%channel_eff)) then
	 allocate(obs(i)%peak%channel_eff(obs(i)%peak%Nc))
	endif 
    obs(i)%peak%channel_eff(:) = obs(i)%table%channel_eff(:) 
	if(.not.allocated(obs(i)%peak%channel_mu)) then
	 allocate(obs(i)%peak%channel_mu(obs(i)%peak%Nc))
	endif 
	if(.not.allocated(obs(i)%peak%channel_w_model)) then
	 allocate(obs(i)%peak%channel_w_model(obs(i)%peak%Nc))
	endif 
	if(.not.allocated(obs(i)%peak%channel_w)) then
	 allocate(obs(i)%peak%channel_w(obs(i)%peak%Nc))
	endif 
	if(.not.allocated(obs(i)%peak%channel_w_corrected_eff)) then
	 allocate(obs(i)%peak%channel_w_corrected_eff(obs(i)%peak%Nc))
	endif 
	if(.not.allocated(obs(i)%peak%channel_systSM)) then
	 allocate(obs(i)%peak%channel_systSM(obs(i)%peak%Nc))
	endif 
	if(.not.allocated(obs(i)%peak%channel_syst)) then
	 allocate(obs(i)%peak%channel_syst(obs(i)%peak%Nc))
	endif 
    do j=1,obs(i)%peak%Nc
     obs(i)%peak%channel_mu(j)=0.0D0
     obs(i)%peak%channel_w(j)=0.0D0
     obs(i)%peak%channel_w_corrected_eff(j)=0.0D0     
     obs(i)%peak%channel_w_model(j)=0.0D0
     obs(i)%peak%channel_systSM(j)=0.0D0 		  
     obs(i)%peak%channel_syst(j)=0.0D0     
    enddo
   endif     
  enddo
   
 end subroutine setup_peak_observables
!--------------------------------------------------------------------
 subroutine read_mutable(t1,skip,col,fullfilename)
!--------------------------------------------------------------------
!--------------------------------------input
  type(mutable) :: t1  
  integer :: skip,col
  character(LEN=*) :: fullfilename
!-----------------------------------internal
  integer :: i,n , file_id_1  
  double precision :: xdummy,xdummy_store
!-------------------------------------------

  file_id_1 = 666
  t1%mu=-1.0D0
  
  open(file_id_1, file=(trim(fullfilename)))
   
  do i=1,skip
   read(file_id_1,*) !skip lines
  enddo 

  xdummy_store = t1%xmin-t1%sep
       
  do i=1,t1%nx
   read(file_id_1,*)xdummy,(t1%mu(i,n),n=1,3)

   ! checks that x are evenly spaced as expected
   if((abs(xdummy-xdummy_store-t1%sep).gt.1.0D-7) &
     &  .or.(abs(xdummy-(t1%xmin+dble(i-1)*t1%sep)).gt.1.0D-7))then
!!    write(*,*)i,t1%id,xdummy,t1%xmin+dble(i-1)*t1%sep
	write(*,*) "Problem with observable ",t1%id
    stop 'error in read_mutable (a1)'
   endif           
   t1%mass(i) = xdummy
!!   write(*,*) xdummy
   xdummy_store=xdummy       
       
  enddo  
   
  if(abs(xdummy-t1%xmax).gt.1.0D-7)stop 'error in read_mutable (a2)'

  close(file_id_1) 
    
 end subroutine read_mutable
!--------------------------------------------------------------------
 subroutine setup_tablelist
!--------------------------------------------------------------------
 implicit none
 integer,allocatable :: tableid(:)
 integer :: i,j,k,n, ntab, firstpeak_index
 logical :: tabused, tab_taken_from_mpredobs

 allocate(tableid(size(obs,dim=1)))
 do i=lbound(tableid,dim=1),ubound(tableid,dim=1)
  tableid(i)=-1
 enddo
 
 ntab=0
 do i=lbound(obs,dim=1),ubound(obs,dim=1)
  tabused=.False.
  do k=lbound(tableid,dim=1),ubound(tableid,dim=1)
   if(obs(i)%table%id.eq.tableid(k)) tabused=.True.
  enddo
  if(.not.tabused) then
   ntab=ntab+1
!--Fill first element which is -1 with the table id   
   do k=lbound(tableid,dim=1),ubound(tableid,dim=1)
    if(tableid(k).eq.-1) then
     tableid(k)=obs(i)%table%id
     exit
    endif 
   enddo
  endif 
 enddo 

 if(allocated(analyses)) deallocate(analyses)
 allocate(analyses(ntab))
 
 do i=lbound(analyses,dim=1),ubound(analyses,dim=1)
  analyses(i)%id = tableid(i)
!--Find observables based on this table
  n=0 
  do j=lbound(obs,dim=1),ubound(obs,dim=1)
   if(analyses(i)%id.eq.obs(j)%table%id) then	 
!----Check if it is a peak observable
    if(obs(j)%obstype.eq.1) then
     n=n+1
    endif
   endif
  enddo
  allocate(analyses(i)%peaks(n))
  analyses(i)%Npeaks=n
  n=0
  firstpeak_index=0 
  tab_taken_from_mpredobs=.False.
  do j=lbound(obs,dim=1),ubound(obs,dim=1)
   if(analyses(i)%id.eq.obs(j)%table%id) then	 
!----Check if it is a peak observable
    if(obs(j)%obstype.eq.1) then
     n=n+1
     analyses(i)%peaks(n) = obs(j)%peak
     analyses(i)%peaks(n)%Higgses = obs(j)%Higgses
	 allocate(analyses(i)%peaks(n)%channel_w_allH( &
&     size(obs(j)%table%channel_w,dim=1),size(obs(j)%table%channel_w,dim=2)))
	 allocate(analyses(i)%peaks(n)%channel_w_corrected_eff_allH( &
&     size(obs(j)%table%channel_w_corrected_eff,dim=1), &
&     size(obs(j)%table%channel_w_corrected_eff,dim=2)))
	 allocate(analyses(i)%peaks(n)%channel_systSM_allH( &
&     size(obs(j)%table%channel_systSM,dim=1),size(obs(j)%table%channel_systSM,dim=2)))
     allocate(analyses(i)%peaks(n)%channel_syst_allH(        &
&     size(obs(j)%table%channel_syst,dim=1),size(obs(j)%table%channel_syst,dim=2)))
     allocate(analyses(i)%peaks(n)%channel_mu_allH(         &
&    size(obs(j)%table%channel_mu,dim=1),size(obs(j)%table%channel_mu,dim=2)))
     analyses(i)%peaks(n)%channel_w_allH = obs(j)%table%channel_w
     analyses(i)%peaks(n)%channel_w_corrected_eff_allH = obs(j)%table%channel_w_corrected_eff
     analyses(i)%peaks(n)%channel_systSM_allH = obs(j)%table%channel_systSM
     analyses(i)%peaks(n)%channel_syst_allH = obs(j)%table%channel_syst
     analyses(i)%peaks(n)%channel_mu_allH = obs(j)%table%channel_mu
     if(n.eq.1) firstpeak_index=j
    elseif(obs(j)%obstype.eq.2) then
     analyses(i)%table = obs(j)%table
     if(.not.allocated(analyses(i)%Higgses)) then
      allocate(analyses(i)%Higgses(size(obs(j)%Higgses,dim=1)))
     endif 
     analyses(i)%Higgses = obs(j)%Higgses
     do k=lbound(analyses(i)%Higgses,dim=1),ubound(analyses(i)%Higgses,dim=1)
      analyses(i)%Higgses(k)%mp_test=1
     enddo 
     tab_taken_from_mpredobs = .True.
    endif 
   endif
  enddo
!-If we did not find a mpred observable, take the table from the first peak observable
 if(.not.tab_taken_from_mpredobs) then
  if(firstpeak_index.gt.0) then
   analyses(i)%table=obs(firstpeak_index)%table
   if(.not.allocated(analyses(i)%Higgses)) then
    allocate(analyses(i)%Higgses(size(obs(firstpeak_index)%Higgses,dim=1)))      
   endif 
   analyses(i)%Higgses=obs(firstpeak_index)%Higgses
  endif    
 endif      
 enddo


 deallocate(tableid)

 end  subroutine setup_tablelist
!--------------------------------------------------------------------
 subroutine setup_observables
!-------------------------------------------------------------------- 
  use usefulbits, only : debug,np,not_a_particle

    implicit none      
  !-----------------------------------internal 
  integer :: i,c, nobs
         
  call initialize_observables
  call setup_peak_observables
  
 !n.b.: setup of mpred observables will be done during the run.
  do i=lbound(obs,dim=1),ubound(obs,dim=1)
   if( count(obs%id.eq.obs(i)%id).ne.1)then
    write(*,*)'the observable id',obs(i)%id,'is repeated'
    stop 'error in setup_observables (c1)'
   endif
  enddo 

  !check to make sure that obs(i)%mutable%particle_x are all particles
  do i=lbound(obs,dim=1),ubound(obs,dim=1)
   if(obs(i)%table%particle_x.eq.not_a_particle)then
    write(*,*)obs(i)%id,'particle_x=not_a_particle.'
    stop 'error in setup_observables (d1)'  
   endif
  enddo

 end subroutine setup_observables
!--------------------------------------------------------------------
 subroutine check_available_Higgses(ii)
!--------------------------------------------------------------------
  implicit none
  integer, intent(in) :: ii
  integer :: i,j
  
  if(.not.allocated(analyses)) then
   stop "Error in subroutine check_available_Higgses: analyses not allocated."
  endif
  
  do i=lbound(analyses(ii)%peaks,dim=1),ubound(analyses(ii)%peaks,dim=1)
   do j=lbound(analyses(ii)%peaks(i)%Higgs_comb,dim=1),    &
&		ubound(analyses(ii)%peaks(i)%Higgs_comb,dim=1)
    if(analyses(ii)%peaks(i)%Higgs_comb(j).ne.0) then
     analyses(ii)%Higgses(analyses(ii)%peaks(i)%Higgs_comb(j))%mp_test=0
    endif 
   enddo
  enddo

 end subroutine check_available_Higgses
!--------------------------------------------------------------------
 subroutine interpolate_mutable(value, mh_in, mutab, steps)
! Returns interpolated values for the signal strength (lower, central, upper value)
! for the given mass mh_in of mutable mutab. "steps" defines the number of interpolation
! steps.
!--------------------------------------------------------------------
  implicit none
   double precision, dimension(3), intent(out) :: value
   double precision, intent(in) :: mh_in
   type(mutable), intent(in) :: mutab
   integer, intent(in) :: steps 
   integer :: i,j, imin, imax
   double precision :: subsep, mh, closestdistance, mh_closest
   double precision, dimension(3) :: mu
      
   if(mh_in.ge.mutab%xmax) then
    value(:) = mutab%mu(ubound(mutab%mu,dim=1),:)
   else if(mh_in.le.mutab%xmin) then 
    value(:) = mutab%mu(lbound(mutab%mu,dim=1),:)   
   else 
    imin = int((mh_in-mutab%xmin)/mutab%sep)+1
    subsep = mutab%sep/steps

    closestdistance = 1000000.0D0
    mh = mutab%mass(imin)
    do j=1, steps
     mh = mh + subsep
     mu(:) = mutab%mu(imin,:) + (mutab%mu(imin+1,:)-mutab%mu(imin,:))*j/steps
     if(abs(mh-mh_in).le.closestdistance) then
      closestdistance = abs(mh-mh_in)
      mh_closest = mh
      value = mu
     endif
    enddo  
    endif
   
 end subroutine interpolate_mutable
!--------------------------------------------------------------------
subroutine deallocate_observables
!--------------------------------------------------------------------
 implicit none
 integer :: x
 
 if(allocated(obs)) then
  do x=lbound(obs,dim=1),ubound(obs,dim=1)
   deallocate(obs(x)%table%mu)
   deallocate(obs(x)%table%mass)  
   deallocate(obs(x)%table%channel_mu)
   deallocate(obs(x)%table%channel_id)  
   deallocate(obs(x)%table%channel_systSM)
   deallocate(obs(x)%table%channel_syst)  
   deallocate(obs(x)%table%channel_w)
   deallocate(obs(x)%table%channel_w_corrected_eff)   
   deallocate(obs(x)%table%channel_eff)
   deallocate(obs(x)%table%channel_eff_ratios)   
   deallocate(obs(x)%table%channel_description)  
   if(allocated(obs(x)%table%Toys_mhobs)) deallocate(obs(x)%table%Toys_mhobs)
   if(allocated(obs(x)%table%Toys_muobs)) deallocate(obs(x)%table%Toys_muobs)
   if(allocated(obs(x)%peak%Higgs_comb)) then
    deallocate(obs(x)%peak%Higgs_comb)
    deallocate(obs(x)%peak%Higgses)
    deallocate(obs(x)%peak%channel_id)
    deallocate(obs(x)%peak%channel_eff)
    deallocate(obs(x)%peak%channel_mu)
    deallocate(obs(x)%peak%channel_w_model)   
    deallocate(obs(x)%peak%channel_w)
    deallocate(obs(x)%peak%channel_w_corrected_eff)    
    deallocate(obs(x)%peak%channel_systSM)
    deallocate(obs(x)%peak%channel_syst)
   endif 
   if(allocated(obs(x)%Higgses)) deallocate(obs(x)%Higgses)    
  enddo
  deallocate(obs)
 endif 
   
end subroutine deallocate_observables
!--------------------------------------------------------------------
end module datatables
!--------------------------------------------------------------------