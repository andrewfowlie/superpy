
!******************************************************
program HBwithSLHA
!
! Run with
!   ./HBwithSLHA npoints <stem>
! where npoints is the number of parameter points you would like to
! look at and each parameter point has a corresponding SLHA file
! e.g. the corresponding SLHA for the 5th point should be found at
!   <stem>.5
!
! Output
! The block HiggsBoundsResults will be added to each SLHA file.
! In addition, the HiggsBounds results will be collected together in 
! the file
!   <stem>-fromHB
!
!******************************************************     
#ifdef NAGf90Fortran
 use F90_UNIX_ENV, only : iargc,getarg
#endif
  implicit none
  integer :: nH,nHplus,HBresult,chan,ncombined
  double precision :: obsratio
  integer :: i,npoints, ios
  integer,parameter :: fileid=78, fileid2=80
  character(len=8) :: istring
  character(len=300) :: inputfilename,outputfilename
  character(len=300) :: stem
  character(LEN=300) :: temp, tmpstring
  integer :: number_args
! used by set_mass_uncertainties
  double precision, allocatable :: dmhneut(:)
  double precision, allocatable :: dmhch(:)
#ifndef NAGf90Fortran
  integer :: iargc
#endif

  nH=3
  nHplus=1

  allocate(dmhneut(nH),dmhch(nHplus))

  number_args = IARGC() 
  
  if( number_args .ne. 2)then
   stop "Incorrect number of arguments given to HBwithSLHA"
  endif

  ! Read arguments into text strings.
  i=1
  temp=""
  call GETARG(i,temp)
  read(temp,*) npoints

  i=i+1  
  temp=""
  call GETARG(i,temp)
  stem = ""
  stem = trim(temp)

  call initialize_HiggsBounds(nH,nHplus,'LandH')

  outputfilename=trim(adjustl(stem))//'-fromHB'

  open(fileid, file=trim(outputfilename))

  do i=1,npoints

   if(i.gt.99999999)stop'need to increase the size of istring in HBwithSLHA'
   write(istring,'(I8)')i

   inputfilename=trim(adjustl(stem))//'.'//trim(adjustl(istring))

 !-Test if input file exists and is non-empty
   open(fileid2, file=inputfilename, form='formatted')
   read(fileid2,'(A)',iostat=ios) tmpstring

   if(ios.eq.0) then
    close(fileid2)   

    call HiggsBounds_input_SLHA(inputfilename)

!    dMhneut = (/2.D0, 0.D0, 0.D0/)
!    dMhch   = (/0.D0/)
!    call HiggsBounds_set_mass_uncertainties(dmhneut,dmhch)

    call run_HiggsBounds( HBresult, chan, obsratio, ncombined )

    !This will add the block HiggsBoundsResults to the SLHA file
    call HiggsBounds_SLHA_output

    !This will collect all the HiggsBounds results together into one file
    write(fileid,*)i,HBresult,chan,obsratio,ncombined

   else
    close(fileid2)     
    call system("rm -f "//inputfilename)    
   endif 


  enddo
  
  close(fileid)

  call finish_HiggsBounds
  deallocate(dmhneut,dmhch)

end program HBwithSLHA
