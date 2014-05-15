!--------------------------------------------------------------------
! This file is part of HiggsSignals (TS 31/01/2013)
!--------------------------------------------------------------------
module combinatorics

 implicit none
 integer, allocatable :: indices(:,:)
 integer, allocatable :: peakindices(:,:,:)

 contains
!--------------------------------------------------------------------
  subroutine create_peakindices(npeaks, nhiggs)
!--------------------------------------------------------------------
 integer, intent(in) :: npeaks, nhiggs
!--Internal parameter:
 integer, allocatable :: Higgs_avail(:,:),Higgs_avail_tmp(:,:)
 integer, allocatable :: peakindices_tmp(:,:,:),peakindices_big(:,:,:) 
 integer :: ncomb, a, b, i, j, N, nphcomb, jj
 logical :: avail, zeros

if(npeaks.le.0) then
 write(*,*) 'WARNING in subroutine create_peakindices:',&
 &   	    'npeaks must be greater than 0!'
else
 call create_indices(nhiggs)
 ncomb=ncombinations(nhiggs)+1

! nphcomb=npeakhiggscombs(nhiggs,npeaks)+1
! allocate(peakindices(nphcomb,npeaks,nhiggs))
! allocate(Higgs_avail(nphcomb,nhiggs))

! We first create very big matrices and reduce them later to the number
! of irreducible combinations: This should be the correct number, otherwise
! there will be a warning later...
 nphcomb = (npeaks+1)**(nhiggs)
 allocate(peakindices_big(nphcomb,npeaks,nhiggs))
 allocate(Higgs_avail(nphcomb,nhiggs))

 
 !Make all Higgses available by setting all entries of Higgs_avail to 1
 do i=lbound(Higgs_avail,dim=1),ubound(Higgs_avail,dim=1)
  do j=lbound(Higgs_avail,dim=2),ubound(Higgs_avail,dim=2)
   Higgs_avail(i,j)=1
  enddo
 enddo
 
 N=1
 do i=1,npeaks
  b=N
  do j=1, ncomb
!  if(i.ne.1) then
 !-Copy all existing matrices to a temporary matrix collection
   allocate(peakindices_tmp(N,npeaks,nhiggs))   
   allocate(Higgs_avail_tmp(N,nhiggs))
   do a=1, N
     call copy_matrices(peakindices_big(a,:,:),peakindices_tmp(a,:,:))
     Higgs_avail_tmp(a,lbound(Higgs_avail_tmp,dim=2):ubound(Higgs_avail_tmp,dim=2)) &
     & =Higgs_avail(a,lbound(Higgs_avail_tmp,dim=2):ubound(Higgs_avail_tmp,dim=2))
   enddo
!  endif     
  
   do a=1, N 
    call check_Higgs_avail(avail, zeros, Higgs_avail_tmp(a,:), indices(j,:))
!	write(*,*) "a=",a,", b=",b," , N=",N
!    write(*,*) "Higgs_avail(a,:) = ",Higgs_avail_tmp(a,:)
!    write(*,*) "indices(j,:) = ",indices(j,:)
!    write(*,*) "avail = ", avail
!    write(*,*) "zeros = ", zeros
    if(avail) then
     if(zeros) then
!      write(*,*) "... filling existing matrix with a row containing zeros..."
!-----Add another row (for this peak) to the existing matrices	     
      peakindices_big(a, i, lbound(peakindices_big,dim=3):ubound(peakindices_big,dim=3)) &
     & = indices(j,lbound(indices,dim=2):ubound(indices,dim=2))
      Higgs_avail(a,lbound(Higgs_avail,dim=2):ubound(Higgs_avail,dim=2)) = &
     & Higgs_avail_tmp(a,lbound(Higgs_avail_tmp,dim=2):ubound(Higgs_avail_tmp,dim=2))
     else
!      write(*,*) "... create new matrix and add a row..."     
      b=b+1
      call copy_matrices(peakindices_tmp(a,:,:),peakindices_big(b,:,:))
      peakindices_big(b,i,lbound(peakindices_big,dim=3):ubound(peakindices_big,dim=3)) &
     & = indices(j,lbound(indices,dim=2):ubound(indices,dim=2))
      Higgs_avail(b,lbound(Higgs_avail,dim=2):ubound(Higgs_avail,dim=2)) = &
     & Higgs_avail_tmp(a,lbound(Higgs_avail_tmp,dim=2):ubound(Higgs_avail_tmp,dim=2))
!DEBUGGING
!      write(*,*) "peakindices(",b,jj,"):"
!      do jj = lbound(peakindices,dim=2),ubound(peakindices,dim=2)
!       write(*,*) peakindices(b,jj,:)
!      enddo 

     endif
    endif
   enddo	
   deallocate(peakindices_tmp)
   deallocate(Higgs_avail_tmp)
  enddo
  N=b

 enddo
 
 deallocate(indices)

!NEW:
 deallocate(Higgs_avail)
 allocate(peakindices(N,npeaks,nhiggs))
 do a=lbound(peakindices,dim=1),ubound(peakindices,dim=1)
  call copy_matrices(peakindices_big(a,:,:),peakindices(a,:,:))
 enddo
 deallocate(peakindices_big) 
 
!- 
 if(N.ne.nphcomb) then
 write(*,*) "WARNING in subroutine create_peakindices:",&
 & " problem with number of possible combinations." 
 write(*,*) "created peakindices for nhiggs/npeak = ",nhiggs,npeaks
 write(*,*) "N = ", N
 write(*,*) "nphcomb = ", nphcomb
 endif
 
 endif
 
 end subroutine create_peakindices
!--------------------------------------------------------------------
 subroutine check_Higgs_avail(avail, zeros, Higgs_avail, Higgs_comb)
!--------------------------------------------------------------------
  logical, intent(out) :: avail, zeros
  integer, dimension(:) :: Higgs_avail
  integer, dimension(:), intent(in) :: Higgs_comb
  integer :: i
 
  !First, check if dimensions match:
  if(size(Higgs_avail).ne.size(Higgs_comb)) then
   write(*,*) 'WARNING in subroutine check_Higgs_avail: dimensions of',&
   &		  'Higgs_avail and Higgs_comb do not match!'
  else
  !-Check, whether the Higgs_comb is still available:
   avail=.True.
   zeros=.True.
   do i=1, size(Higgs_comb)
    if(avail) then
     if(Higgs_comb(i).ne.0) then
      zeros=.False.
      if(Higgs_avail(Higgs_comb(i)).eq.0) avail=.False.
     endif
    endif
   enddo
  
   if(avail.and.(.not.zeros)) then
   !-This combination is available. Now make the used Higgses unavailable
   !-for other combinations
    do i=1, size(Higgs_comb)
     if(Higgs_comb(i).ne.0) Higgs_avail(Higgs_comb(i)) = 0
    enddo       
   endif
  endif 
  
 end subroutine check_Higgs_avail
!--------------------------------------------------------------------
 subroutine copy_matrices(matrix_in, matrix_out)
!--------------------------------------------------------------------
 integer, dimension(:,:), intent(in) :: matrix_in
 integer, dimension(:,:), intent(out) :: matrix_out
 integer :: a,b
 
 
   if(size(matrix_in,dim=1).ne.size(matrix_out,dim=1).or. &
&     size(matrix_in,dim=2).ne.size(matrix_out,dim=2)) then
   write(*,*) 'WARNING in subroutine copy matrices: dimensions of',&
   &		  'matrices do not match!'
   else
    do a=lbound(matrix_in,dim=1),ubound(matrix_in,dim=1)
     do b=lbound(matrix_in,dim=2),ubound(matrix_in,dim=2)
	  matrix_out(a,b)=matrix_in(a,b)
     enddo
    enddo    
   endif 
 
 end subroutine copy_matrices
!--------------------------------------------------------------------
  subroutine create_indices(nhiggs)
!--------------------------------------------------------------------
  integer, intent(in) :: nhiggs
  integer :: ncombs, a, i, j, k, x, ii, m, leader
  integer, allocatable :: s(:)
  
  ! Calculate the number of possible (irreducible) Higgs combinations
  ncombs=ncombinations(nhiggs)+1
  
!  write(*,*) 'ncombs = ', ncombs
  
  allocate(indices(ncombs,nhiggs))
  ! Fill the matrix indices with all possible Higgs combinations
  ! Fill the first row of indices with zeros
  do i=lbound(indices,dim=2),ubound(indices,dim=2)
   indices(1,i)=0
  enddo
  
    
! Fill the matrix column by column (i,j)
  do j=1,nhiggs
   a=1
! First column
   if(j.eq.1) then
    do i=1,nhiggs
     do ii=1,2**(i-1)
      a=a+1
      indices(a,j)=i
     enddo
    enddo 
! Preceding columns
   else
    a=1
	do while (a < ncombs)
     a=a+1
     leader = indices(a,j-1)
     indices(a,j)=0
     if(leader.ne.1) then
      do i=1,leader-1
       do ii=1,2**(i-1)
        a=a+1
        indices(a,j)=i
       enddo
      enddo  
     endif 
    enddo
   endif 
  enddo
  
!  write(*,*) "Indices:"
!  do i=1,ncombs
!   write(*,*) indices(i,:)
!  enddo 
!  if(a.ne.ncombs) write(*,*) "WARNING in subroutine create indices:",&
!  				" Higgs combinations do not match! a=",a,", ncomb=",ncombs

  end subroutine create_indices
!--------------------------------------------------------------------
  function ncombinations(nhiggs)
!--------------------------------------------------------------------
  integer, intent(in) :: nhiggs
  integer :: ncombinations, x
  integer, allocatable :: s(:)
  integer :: kk
  ! Calculate the number of possible (irreducible) Higgs combinations
  allocate(s(nhiggs+1))  
  s(1)=0
  do x=1,nhiggs	
   s(x+1)=s(x)+(x-1)
  enddo
  ncombinations=sum(s)+nhiggs+1 
!NEW
  deallocate(s)
  
  ncombinations=0
  do kk=0,nhiggs-1
   ncombinations = ncombinations + 2**kk 
  enddo

  end function ncombinations
!--------------------------------------------------------------------
!  function nsubcombinations(nhiggs, ncomb)
!!--------------------------------------------------------------------
!  integer, intent(in) :: nhiggs, ncomb
!  integer :: nsubcombinations
!  integer :: i
!  
!  if((nhiggs.lt.ncomb).or.(nhiggs.le.0).or.(ncomb.lt.0)) then
!   write(*,*) "WARNING in function nsubcombinations: nhiggs or ncomb not valid!"
!  else
!   if(ncomb.eq.0) then
!    nsubcombinations = 1
!   else if(ncomb.eq.1) then
!    nsubcombinations = nhiggs
!   else
!    nsubcombinations=0
!	do i=ncomb,nhiggs
!	 nsubcombinations=nsubcombinations+(i-ncomb+1)          
!	enddo 
!   endif 
!  endif
!  end function nsubcombinations  
!--------------------------------------------------------------------
!  recursive function npeakhiggscombs(nhiggs, npeaks) result(N)
!!--------------------------------------------------------------------
!  integer, intent(in) :: nhiggs, npeaks
!  integer :: j,k,N
!
!!  write(*,*) "Start npeakhiggscombs with nhiggs=",nhiggs,", npeaks=",npeaks
!
!  if(npeaks.lt.1) then
!   write(*,*) "WARNING in function npeakhiggscombs: npeak not valid!"
!   N=0
!  else if(nhiggs.lt.1) then
!   N=1
!  else 
!   if(npeaks.eq.1) then
!    N=ncombinations(nhiggs)
!   else if(nhiggs.eq.1) then
!    N=npeaks+1
!   else
!    N=0 
!    do j=0, nhiggs
!     do k=1,nsubcombinations(nhiggs,j)
!      N=N+npeakhiggscombs(nhiggs-j,npeaks-1)
!     enddo 
!    enddo 
!   endif
!  endif
!
!  end function npeakhiggscombs
!--------------------------------------------------------------------
end module combinatorics
!--------------------------------------------------------------------