!********************************************************
module SLHA_manip
!********************************************************
 use string_manip

 implicit none
 private 
 public :: readSLHAfile,finishwithSLHA, &
    &      get_mass,get_totaldecaywidth,get_twobodybranchingratio, &
    &      get_SPhenocrosssection,get_SPhenocrosssectionCMenergy,&
    &      get_threebodybranchingratio,get_modsel,get_HBresults_obsratio,&
    &      get_HiggsCouplingsBosons,get_HiggsCouplingsFermions, &
    &      writeSLHAfile_except,check_validity,get_HBresults_channel_id, &
    &      get_mass_uncertainty


 type array_of_strings      
  character(len=200) :: line
  integer :: id       
 end type

 type array_of_integers_and_dble     
  double precision :: val
  integer :: id       
 end type

 type(array_of_strings), allocatable :: contents_of_input_file(:)
 type(array_of_strings), allocatable :: edited_contents_of_input_file(:)
 integer :: n_input_lines  

 type(array_of_strings), allocatable :: names_of_blocks_and_decays(:)
 type(array_of_integers_and_dble), allocatable :: store_masses(:)
 type(array_of_integers_and_dble), allocatable :: store_mass_uncertainties(:) 

 contains

 
 !********************************************************
 subroutine readSLHAfile(fileid)
 !reads in contents of SLHA file, and stores it in 'contents_of_input_file'
 ! strips off comments and tabs and
 ! stores it in 'edited_contents_of_input_file' line-by-line 
 !********************************************************
  implicit none 
  integer,intent(in) :: fileid
  integer :: i,a
  integer :: stringlength
  character(len=200) :: temp_line
  
  !finds number of lines to read in
  n_input_lines=getSLHAfilelength(fileid)
  allocate(contents_of_input_file(n_input_lines)) 
  allocate(edited_contents_of_input_file(n_input_lines))   

  !reads each line in to 'contents_of_input_file'
  do i=1,n_input_lines
   read(fileid,'(a)')contents_of_input_file(i)%line
  enddo

  !this gets rid of the comments and tabs and puts the result into 'edited_contents_of_input_file'
  do i=1,n_input_lines
   temp_line=contents_of_input_file(i)%line
   stringlength=len(temp_line)
   temp_line=strip_off_comment(temp_line,stringlength)  

   !-------------
   ! should not be needed, since SLHA files should not contain tabs
   ! but it's safer to check
   a=1
   do while (a.gt.0)  !replace all tab characters with spaces
    a=INDEX(temp_line,'	')
    temp_line=temp_line(:a-1)//'   '//temp_line(a+1:)
   enddo
   !-------------

   edited_contents_of_input_file(i)%line=temp_line
  enddo 

  call store_line_numbers_of_blocks_and_decays
  call fill_store_masses
  call fill_store_mass_uncertainties

  !do i=1,ubound(names_of_blocks_and_decays,dim=1)
  ! write(*,*)'hello block',names_of_blocks_and_decays(i)%id, &
  !  & trim(adjustl(names_of_blocks_and_decays(i)%line))
  !enddo
  !stop'here for now (subroutine readSLHAfile)'

 end subroutine readSLHAfile
 !********************************************************
 subroutine store_line_numbers_of_blocks_and_decays
 !********************************************************
  integer :: x,n,xtot
  character(len=200) :: temp_line
  integer :: stringlength,stringlength2
  character(len=50),allocatable :: col(:)
  integer :: number_of_blocks_and_decays
  integer, allocatable :: blockordecaystart(:)
  character(len=1) :: firstchar

  allocate(col(2))
  xtot=ubound(edited_contents_of_input_file,dim=1)
  allocate(blockordecaystart(xtot))
  stringlength =len(edited_contents_of_input_file%line)
  stringlength2=len(col)

  blockordecaystart=0

  do x=1,xtot
   
    temp_line=edited_contents_of_input_file(x)%line
    temp_line=adjustl(temp_line)
    if(trim(temp_line).ne.'')then
     firstchar=temp_line(1:1)
     if(     (firstchar.eq.'D').or.(firstchar.eq.'d')  & 
      &  .or.(firstchar.eq.'b').or.(firstchar.eq.'B')  )then
          call split_into_col(temp_line,stringlength,col,stringlength2)

          col(1) = strtolcase( col(1) ,stringlength2)
  
          if((col(1).eq.'block').or.(col(1).eq.'decay'))then
            blockordecaystart(x)=1
          endif
      endif
     endif
  enddo

  number_of_blocks_and_decays=sum(blockordecaystart)

  allocate(names_of_blocks_and_decays(number_of_blocks_and_decays))

  n=0
  do x=1,xtot
   if(blockordecaystart(x).eq.1)then
    temp_line=edited_contents_of_input_file(x)%line

    call split_into_col(temp_line,stringlength,col,stringlength2)

    n=n+1
    names_of_blocks_and_decays(n)%id=x   
    names_of_blocks_and_decays(n)%line=strtolcase(col(2),stringlength2)
   endif
  enddo

  deallocate(col)
  deallocate(blockordecaystart)
 end subroutine store_line_numbers_of_blocks_and_decays
                                                                                                             
 !********************************************************
 subroutine writeSLHAfile_except(fileid,name_of_block1,name_of_block2,name_of_block3)
 ! writes the inputted data from an SLHA file to the file with
 ! the filehandle=fileid. Includes original comments and tabs
 ! does not write out the block with name 'name_of_block'
 !********************************************************
  implicit none 
  integer,intent(in) :: fileid
  character(len=*),intent(in),optional ::name_of_block1,name_of_block2,name_of_block3
  integer :: blocklines1(2),blocklines2(2),blocklines3(2),blocklines_temp(2)
  integer :: i
  
  if(present(name_of_block1))then
   blocklines1=line_numbers_of_block(name_of_block1)
  else
   blocklines1=0
  endif
  if(present(name_of_block2))then
   blocklines2=line_numbers_of_block(name_of_block2)
  else
   blocklines2=0
  endif
  if(present(name_of_block3))then
   blocklines3=line_numbers_of_block(name_of_block3)
  else
   blocklines3=0
  endif
  

  if(    blocklines1(1).eq.blocklines2(1))then ! they are the same block
    blocklines1=0
  elseif(blocklines1(1).gt.blocklines2(1))then ! swap them
    blocklines_temp=blocklines1
    blocklines1=blocklines2
    blocklines2=blocklines_temp
  endif

  if(    blocklines1(1).eq.blocklines3(1))then ! they are the same block
    blocklines1=0
  elseif(blocklines1(1).gt.blocklines3(1))then ! swap them
    blocklines_temp=blocklines1
    blocklines1=blocklines3
    blocklines3=blocklines_temp
  endif

  if(    blocklines2(1).eq.blocklines3(1))then ! they are the same block
    blocklines2=0
  elseif(blocklines2(1).gt.blocklines3(1))then ! swap them
    blocklines_temp=blocklines2
    blocklines2=blocklines3
    blocklines3=blocklines_temp
  endif


  do i=1,ubound(contents_of_input_file,dim=1)

    if(  (                  i.lt.blocklines1(1)                ) .or. &
      &  ( ( i.gt.blocklines1(2) ).and.( i.lt.blocklines2(1) ) ) .or. &
      &  ( ( i.gt.blocklines2(2) ).and.( i.lt.blocklines3(1) ) ) .or. &
      &  (                  i.gt.blocklines3(2)                ) )then
     write(fileid,'(a)')trim(contents_of_input_file(i)%line)
    endif

  enddo
 
 end subroutine writeSLHAfile_except
 !********************************************************
 subroutine finishwithSLHA
 !********************************************************
  implicit none 
   if(allocated(contents_of_input_file)) deallocate(contents_of_input_file)
   if(allocated(edited_contents_of_input_file)) deallocate(edited_contents_of_input_file)
   if(allocated(names_of_blocks_and_decays)) deallocate(names_of_blocks_and_decays)
   if(allocated(store_masses)) deallocate(store_masses)
   if(allocated(store_mass_uncertainties)) deallocate(store_mass_uncertainties)   
 end subroutine finishwithSLHA

 !********************************************************
 function get_totaldecaywidth(particlePDGcode)
 !********************************************************
  implicit none
  integer,intent(in) :: particlePDGcode
  integer :: particlePDGcode_fromblock
  integer :: blocklines(2)
  double precision :: get_totaldecaywidth
  character(len=50),allocatable :: col(:)
  character(len=200) :: temp_line
  integer :: stringlength,stringlength2
  integer :: i

  allocate(col(3))
  get_totaldecaywidth=0.0D0

  blocklines=line_numbers_of_decay(particlePDGcode)

  if(blocklines(1).le.0)then
    write(*,*)'Warning: was not able to read the total decay width of particle with PDG code:',particlePDGcode
  else

    i=blocklines(1)
    temp_line=edited_contents_of_input_file(i)%line

       stringlength =len(temp_line)
       stringlength2=len(col) 

       call split_into_col(temp_line,stringlength,col,stringlength2)

       call saferead_int(col(2),particlePDGcode_fromblock)
       if(particlePDGcode.ne.particlePDGcode_fromblock)then
         stop'problem in function get_particledecaywidth(2)'
       endif
  
       call saferead_dble(col(3),get_totaldecaywidth)
  endif

  deallocate(col)

 end function get_totaldecaywidth
 !********************************************************
 function get_SPhenocrosssectionCMenergy(particlePDGcode_1,particlePDGcode_2,pol_1,pol_2,incISR)
 !get the centre of mass of the SPheno cross section
 !********************************************************
  implicit none
  integer,intent(in) :: particlePDGcode_1,particlePDGcode_2,incISR
  double precision,intent(in) :: pol_1,pol_2
  integer :: incISR_fromblock
  double precision :: pol_1_fromblock,pol_2_fromblock
  integer :: blocklines(2),particlePDGcode_req(2),particlePDGcode_fromblock(2)
  double precision :: get_SPhenocrosssectionCMenergy
  character(len=50),allocatable :: col(:)
  character(len=200) :: temp_line
  integer :: stringlength,stringlength2
  integer :: i

  allocate(col(7))
  get_SPhenocrosssectionCMenergy=0.0D0
  particlePDGcode_req(1)=particlePDGcode_1
  particlePDGcode_req(2)=particlePDGcode_2

  blocklines=line_numbers_of_block('SPhenoCrossSections')

  if(minval(blocklines).le.0)then
    write(*,*)'Warning: was not able to read the block SPhenoCrossSections'
 
  else
    i=blocklines(1)+1
    temp_line=edited_contents_of_input_file(i)%line

    stringlength =len(temp_line)
    stringlength2=len(col) 

    call split_into_col(temp_line,stringlength,col,stringlength2)

    call saferead_int(col(2),particlePDGcode_fromblock(1))
    call saferead_int(col(3),particlePDGcode_fromblock(2))
    call saferead_dble(col(5),pol_1_fromblock)
    call saferead_dble(col(6),pol_2_fromblock)
    call saferead_int(col(7),incISR_fromblock)

    if(trim(col(1))           .ne. 'XS'                       )then
      stop'problem in function get_SPhenocrosssectionCMenergy(1)'
    elseif( .not.same_particles(particlePDGcode_req,particlePDGcode_fromblock) )then
      stop'problem in function get_SPhenocrosssectionCMenergy(2)'
    elseif( pol_1             .ne. pol_1_fromblock            )then
      stop'problem in function get_SPhenocrosssectionCMenergy(4)'
    elseif( pol_2             .ne. pol_2_fromblock            )then
      stop'problem in function get_SPhenocrosssectionCMenergy(5)'
    elseif( incISR_fromblock  .ne. incISR_fromblock           )then
      stop'problem in function get_SPhenocrosssectionCMenergy(6)'
    endif

    call saferead_dble(col(4),get_SPhenocrosssectionCMenergy)
  endif

  deallocate(col)

 end function get_SPhenocrosssectionCMenergy
 !********************************************************
 function get_twobodybranchingratio(particlePDGcode_parent,particlePDGcode_da1,particlePDGcode_da2)
 !********************************************************
  implicit none
  integer,intent(in) :: particlePDGcode_parent,particlePDGcode_da1,particlePDGcode_da2
  integer :: particlePDGcode_fromblock_da(2),particlePDGcode_req(2)
  integer :: blocklines(2)
  double precision :: get_twobodybranchingratio
  character(len=50),allocatable :: col(:)
  character(len=200) :: temp_line
  !integer :: stringlength,stringlength2
  integer :: i
  integer :: nda,ios
  double precision :: br_fromblock

  allocate(col(4))
  get_twobodybranchingratio=0.0D0
  blocklines=line_numbers_of_decay(particlePDGcode_parent)
  particlePDGcode_req(1)=particlePDGcode_da1
  particlePDGcode_req(2)=particlePDGcode_da2

  if(minval(blocklines).le.0)then
    write(*,*)'Warning: was not able to read the decay widths of particle with PDG code',particlePDGcode_parent
  else

    do i=1+blocklines(1),blocklines(2)
     temp_line=edited_contents_of_input_file(i)%line
     temp_line=adjustl(temp_line)
 
     if(trim(temp_line).ne.'')then
       read(temp_line,*,iostat=ios)br_fromblock,nda, &
                      & particlePDGcode_fromblock_da(1),particlePDGcode_fromblock_da(2)
       if(ios.eq.0)then
        if(nda.eq.2)then
         if(same_particles(particlePDGcode_req,particlePDGcode_fromblock_da))then
           get_twobodybranchingratio=br_fromblock
         endif
        endif
       endif

       !old method: very safe but too slow
       !stringlength =len(temp_line)
       !stringlength2=len(col) 

       !call split_into_col(temp_line,stringlength,col,stringlength2)

       !call saferead_int(col(2),nda)
       !if(nda.eq.2)then
       ! call saferead_int(col(3),particlePDGcode_fromblock_da(1))
       ! call saferead_int(col(4),particlePDGcode_fromblock_da(2))

       ! if(same_particles(particlePDGcode_req,particlePDGcode_fromblock_da))then

       !   call saferead_dble(col(1),get_twobodybranchingratio)
    
       ! endif
       !endif
     endif
    enddo
  endif

  deallocate(col)

 end function get_twobodybranchingratio
 !********************************************************
 function get_threebodybranchingratio(particlePDGcode_parent,particlePDGcode_da1,particlePDGcode_da2,particlePDGcode_da3)
 !********************************************************
  implicit none
  integer,intent(in) :: particlePDGcode_parent,particlePDGcode_da1,particlePDGcode_da2,particlePDGcode_da3
  integer :: particlePDGcode_fromblock_da(3),particlePDGcode_req(3)
  integer :: blocklines(2)
  double precision :: get_threebodybranchingratio
  character(len=50),allocatable :: col(:)
  character(len=200) :: temp_line
  integer :: stringlength,stringlength2
  integer :: i
  integer :: nda

  allocate(col(5))
  get_threebodybranchingratio=0.0D0
  blocklines=line_numbers_of_decay(particlePDGcode_parent)
  particlePDGcode_req(1)=particlePDGcode_da1
  particlePDGcode_req(2)=particlePDGcode_da2
  particlePDGcode_req(3)=particlePDGcode_da3

  if(minval(blocklines).le.0)then
    write(*,*)'Warning: was not able to read the decay widths of particle with PDG code',particlePDGcode_parent
  else

    do i=1+blocklines(1),blocklines(2)
     temp_line=edited_contents_of_input_file(i)%line

     if(trim(temp_line).ne.'')then

       stringlength =len(temp_line)
       stringlength2=len(col) 

       call split_into_col(temp_line,stringlength,col,stringlength2)

       call saferead_int(col(2),nda)
       if(nda.eq.3)then
        call saferead_int(col(3),particlePDGcode_fromblock_da(1))
        call saferead_int(col(4),particlePDGcode_fromblock_da(2))
        call saferead_int(col(5),particlePDGcode_fromblock_da(3))

        if(same_particles(particlePDGcode_req,particlePDGcode_fromblock_da))then

          call saferead_dble(col(1),get_threebodybranchingratio)
    
        endif
       endif
     endif
    enddo
  endif

  deallocate(col)

 end function get_threebodybranchingratio

 !********************************************************
 function get_SPhenocrosssection(particlePDGcode_da1,particlePDGcode_da2)
 !********************************************************
  implicit none
  integer,intent(in) :: particlePDGcode_da1,particlePDGcode_da2
  integer :: particlePDGcode_fromblock_da(2),particlePDGcode_req(2)
  integer :: blocklines(2)
  double precision :: get_SPhenocrosssection
  character(len=50),allocatable :: col(:)
  character(len=200) :: temp_line
  integer :: stringlength,stringlength2
  integer :: i
  integer :: nda

  allocate(col(4))
  get_SPhenocrosssection=0.0D0
  blocklines=line_numbers_of_block('SPhenoCrossSections')
  particlePDGcode_req(1)=particlePDGcode_da1
  particlePDGcode_req(2)=particlePDGcode_da2

  if(minval(blocklines).le.0)then
    write(*,*)'Warning: was not able to read the block SPhenoCrossSections'
  else

    do i=1+blocklines(1),blocklines(2)
     temp_line=edited_contents_of_input_file(i)%line

     if(trim(temp_line).ne.'')then

       stringlength =len(temp_line)
       stringlength2=len(col) 

       call split_into_col(temp_line,stringlength,col,stringlength2)

       call saferead_int(col(2),nda)
       if(nda.eq.2)then
         call saferead_int(col(3),particlePDGcode_fromblock_da(1))
         call saferead_int(col(4),particlePDGcode_fromblock_da(2))

         if(same_particles(particlePDGcode_req,particlePDGcode_fromblock_da))then

           call saferead_dble(col(1),get_SPhenocrosssection)
    
         endif
       endif
     endif
    enddo
  endif

  deallocate(col)

 end function get_SPhenocrosssection

 !********************************************************
 function get_HiggsCouplingsBosons(particlePDGcode_1,particlePDGcode_2, &
     &   particlePDGcode_3,particlePDGcode_4)
 !********************************************************
 !doesn't matter if particles or antiparticles are used: takes abs() of the PDG number anyway
  implicit none
  integer,intent(in) :: particlePDGcode_1,particlePDGcode_2,particlePDGcode_3
  integer,optional,intent(in) :: particlePDGcode_4
  integer,allocatable :: particlePDGcode_req(:)
  integer,allocatable :: particlePDGcode_fromblock(:)
  integer :: blocklines(2)
  double precision :: get_HiggsCouplingsBosons
  character(len=50),allocatable :: col(:)
  character(len=200) :: temp_line
  !integer :: stringlength,stringlength2
  integer :: i,n
  integer :: np,np_req
  integer :: ios
  double precision :: coup_fromblock

  if(present(particlePDGcode_4))then
    np_req=4
  else
    np_req=3
  endif
 
  allocate(col(np_req+2))
  allocate(particlePDGcode_req(np_req))
  allocate(particlePDGcode_fromblock(np_req))

  get_HiggsCouplingsBosons=0.0D0
  blocklines=line_numbers_of_block('HiggsBoundsInputHiggsCouplingsBosons')
  particlePDGcode_req(1)=particlePDGcode_1
  particlePDGcode_req(2)=particlePDGcode_2
  particlePDGcode_req(3)=particlePDGcode_3

  if(np_req.eq.4)then
    particlePDGcode_req(4)=particlePDGcode_4
  endif

  particlePDGcode_req=abs(particlePDGcode_req)

  if(minval(blocklines).le.0)then
   write(*,*)'Problem reading block HiggsBoundsInputHiggsCouplingsBosons:'
   write(*,*)'Note that this is a non-standard SLHA block and may'
   write(*,*)'need to be entered into the SLHA input file by hand.'
  else

    do i=1+blocklines(1),blocklines(2)
     temp_line=edited_contents_of_input_file(i)%line
     temp_line=adjustl(temp_line)

     if(trim(temp_line).ne.'')then

       read(temp_line,*,iostat=ios)coup_fromblock,np, &
                      & (particlePDGcode_fromblock(n),n=1,np_req)
       particlePDGcode_fromblock=abs(particlePDGcode_fromblock)
       if(ios.eq.0)then
         if(np.eq.np_req)then
           if(same_particles(particlePDGcode_req,particlePDGcode_fromblock))then
             get_HiggsCouplingsBosons=coup_fromblock
           endif
         endif
       endif

       !old method: safe but takes too long
       !stringlength =len(temp_line)
       !stringlength2=len(col) 

       !call split_into_col(temp_line,stringlength,col,stringlength2)
  
       !call saferead_int(col(2),np)

       !if(np.eq.np_req)then

       ! call saferead_int(col(3),particlePDGcode_fromblock(1))
       ! call saferead_int(col(4),particlePDGcode_fromblock(2))
       ! call saferead_int(col(5),particlePDGcode_fromblock(3))
       ! if(np_req.eq.4)then
       !   call saferead_int(col(6),particlePDGcode_fromblock(4))
       ! endif
       ! particlePDGcode_fromblock=abs(particlePDGcode_fromblock)

       ! if(same_particles(particlePDGcode_req,particlePDGcode_fromblock))then

       !   call saferead_dble(col(1),get_HiggsCouplingsBosons)
    
       ! endif

       !endif
     endif
    enddo
  endif


  deallocate(col)
  deallocate(particlePDGcode_req)
  deallocate(particlePDGcode_fromblock)

 end function get_HiggsCouplingsBosons

 !********************************************************
 function get_HiggsCouplingsFermions(particlePDGcode_1,particlePDGcode_2,particlePDGcode_3)
 !********************************************************
 !doesn't matter if particles or antiparticles are used: takes abs() of the PDG number anyway
  implicit none
  integer,intent(in) :: particlePDGcode_1,particlePDGcode_2,particlePDGcode_3
  integer :: particlePDGcode_req(3)
  integer :: particlePDGcode_fromblock(3)
  integer :: blocklines(2)
  double precision :: get_HiggsCouplingsFermions(2)
  character(len=50),allocatable :: col(:)
  character(len=200) :: temp_line
  !integer :: stringlength,stringlength2
  integer :: i
  integer :: np
  integer :: ios
  double precision :: coup_fromblock(2)

  allocate(col(6))
  get_HiggsCouplingsFermions=0.0D0
  blocklines=line_numbers_of_block('HiggsBoundsInputHiggsCouplingsFermions')
  particlePDGcode_req(1)=particlePDGcode_1
  particlePDGcode_req(2)=particlePDGcode_2
  particlePDGcode_req(3)=particlePDGcode_3
  particlePDGcode_req=abs(particlePDGcode_req)

  if(minval(blocklines).le.0)then
   write(*,*)'Problem reading block HiggsBoundsInputHiggsCouplingsFermions:'
   write(*,*)'Note that this is a non-standard SLHA block and may'
   write(*,*)'need to be entered into the SLHA input file by hand.'
  else

    do i=1+blocklines(1),blocklines(2)
     temp_line=edited_contents_of_input_file(i)%line
     temp_line=adjustl(temp_line)

     if(trim(temp_line).ne.'')then

       read(temp_line,*,iostat=ios)coup_fromblock(1),coup_fromblock(2),np, &
                      & particlePDGcode_fromblock(1),particlePDGcode_fromblock(2), &
                      & particlePDGcode_fromblock(3)
       particlePDGcode_fromblock=abs(particlePDGcode_fromblock)
       if(ios.eq.0)then
        if(np.eq.3)then
         if(same_particles(particlePDGcode_req,particlePDGcode_fromblock))then
           get_HiggsCouplingsFermions=coup_fromblock
         endif
        endif
       endif
       !old method: safe but takes too long
       !stringlength =len(temp_line)
       !stringlength2=len(col) 

       !call split_into_col(temp_line,stringlength,col,stringlength2)

       !call saferead_int(col(3),np)
       !if(np.eq.3)then
       ! call saferead_int(col(4),particlePDGcode_fromblock(1))
       ! call saferead_int(col(5),particlePDGcode_fromblock(2))
       ! call saferead_int(col(6),particlePDGcode_fromblock(3))
       ! particlePDGcode_fromblock=abs(particlePDGcode_fromblock)

       ! if(same_particles(particlePDGcode_req,particlePDGcode_fromblock))then

       !   call saferead_dble(col(1),get_HiggsCouplingsFermions(1))
       !   call saferead_dble(col(2),get_HiggsCouplingsFermions(2))   
       ! endif
       !endif
     endif
    enddo
  endif


  deallocate(col)

 end function get_HiggsCouplingsFermions

 !********************************************************
 subroutine fill_store_masses
 !********************************************************
  implicit none
  integer :: blocklines(2)
  integer :: i,n
  integer :: stringlength,stringlength2
  character(len=200) :: temp_line
  character(len=50),allocatable :: col(:)

  allocate(col(2))

  blocklines=line_numbers_of_block('mass')

  if(minval(blocklines).le.0)then
   write(*,*)'Problem reading the block MASS'
   write(*,*)'This block is necessary input to HiggsBounds.' 
   allocate(store_masses(1))
   store_masses(1)%id =0
   store_masses(1)%val=0.0d0
  else
   n=0
   do i=1+blocklines(1),blocklines(2)
    temp_line=edited_contents_of_input_file(i)%line
 
    if(trim(temp_line).ne.'')then !count the lines with something in them
      n=n+1
    endif

   enddo

   allocate(store_masses(n))

   n=0
   do i=1+blocklines(1),blocklines(2)
    temp_line=edited_contents_of_input_file(i)%line

    if(trim(temp_line).ne.'')then
     n=n+1
     !-----------------------
     stringlength =len(temp_line)
     stringlength2=len(col) 

     call split_into_col(temp_line,stringlength,col,stringlength2)

     call saferead_int( col(1),store_masses(n)%id)
     call saferead_dble(col(2),store_masses(n)%val)

    endif
   enddo

   deallocate(col)
  endif

 end subroutine fill_store_masses
 !********************************************************
 subroutine fill_store_mass_uncertainties
 !********************************************************
  implicit none
  integer :: blocklines(2)
  integer :: i,n
  integer :: stringlength,stringlength2
  character(len=200) :: temp_line
  character(len=50),allocatable :: col(:)

  allocate(col(2))

  blocklines=line_numbers_of_block('dmass')

  if(minval(blocklines).le.0)then
   write(*,*)'Problem reading the block DMASS'
   write(*,*)'This block is optional input to HiggsBounds.' 
   write(*,*)'Mass uncertainties will be set to zero.'
   allocate(store_mass_uncertainties(1))
   store_mass_uncertainties(1)%id =0
   store_mass_uncertainties(1)%val=0.0d0
  else
   n=0
   do i=1+blocklines(1),blocklines(2)
    temp_line=edited_contents_of_input_file(i)%line
    if(trim(temp_line).ne.'')then !count the lines with something in them
      n=n+1
    endif

   enddo

   allocate(store_mass_uncertainties(n))

   n=0
   do i=1+blocklines(1),blocklines(2)
    temp_line=edited_contents_of_input_file(i)%line

    if(trim(temp_line).ne.'')then
     n=n+1
     !-----------------------
     stringlength =len(temp_line)
     stringlength2=len(col) 

     call split_into_col(temp_line,stringlength,col,stringlength2)

     call saferead_int( col(1),store_mass_uncertainties(n)%id)
     call saferead_dble(col(2),store_mass_uncertainties(n)%val)

    endif
   enddo

   deallocate(col)
  endif

 end subroutine fill_store_mass_uncertainties 
 !********************************************************
 function get_mass(particlePDGcode)
 !********************************************************
  implicit none
  integer,intent(in) :: particlePDGcode
  integer :: i
  double precision :: get_mass

  get_mass=-1.0D0

  do i=1,ubound(store_masses,dim=1)
     if(particlePDGcode.eq.store_masses(i)%id)then
        get_mass=store_masses(i)%val
     endif
  enddo

 end function get_mass
 !********************************************************
 function get_mass_uncertainty(particlePDGcode)
 !********************************************************
  implicit none
  integer,intent(in) :: particlePDGcode
  integer :: i
  double precision :: get_mass_uncertainty

  get_mass_uncertainty=0.0D0

  do i=1,ubound(store_mass_uncertainties,dim=1)
     if(particlePDGcode.eq.store_mass_uncertainties(i)%id)then
        get_mass_uncertainty=store_mass_uncertainties(i)%val
     endif
  enddo

 end function get_mass_uncertainty
 !********************************************************
 function get_modsel(modsel_id)
 !********************************************************
  implicit none
  integer,intent(in) :: modsel_id
  integer :: modsel_id_fromblock
  integer :: blocklines(2)
  integer :: i
  integer :: stringlength,stringlength2
  character(len=200) :: temp_line
  integer :: get_modsel
  character(len=50),allocatable :: col(:)

  allocate(col(2))
  get_modsel=0

  blocklines=line_numbers_of_block('modsel')
  if(minval(blocklines).le.0)then
   write(*,*)'Problem reading the block modsel'
  else

    do i=1+blocklines(1),blocklines(2)
     temp_line=edited_contents_of_input_file(i)%line

     if(trim(temp_line).ne.'')then

      stringlength =len(temp_line)
      stringlength2=len(col) 

      call split_into_col(temp_line,stringlength,col,stringlength2)

      call saferead_int(col(1),modsel_id_fromblock)

      if(modsel_id.eq.modsel_id_fromblock)then
       call saferead_int(col(2),get_modsel)
      endif
     endif
    enddo
  endif

  deallocate(col)

 end function get_modsel 
 !********************************************************
 function get_HBresults_obsratio(channel_rank)
 !********************************************************
  implicit none
  integer :: id_fromblock(2)
  integer :: blocklines(2)
  integer :: i
  integer, intent(in) :: channel_rank
  integer :: stringlength,stringlength2
  character(len=200) :: temp_line
  double precision :: get_HBresults_obsratio
  character(len=50),allocatable :: col(:)

  allocate(col(3))
  get_HBresults_obsratio=-1.0D0

  blocklines=line_numbers_of_block('HiggsBoundsResults')
  if(minval(blocklines).le.0)then
   write(*,*)'Problem reading the block HiggsBoundsResults'
  else

    do i=1+blocklines(1),blocklines(2)
     temp_line=edited_contents_of_input_file(i)%line

     if(trim(temp_line).ne.'')then

      stringlength =len(temp_line)
      stringlength2=len(col) 

      call split_into_col(temp_line,stringlength,col,stringlength2)

      call saferead_int(col(1),id_fromblock(1))
      call saferead_int(col(2),id_fromblock(2))
      if((id_fromblock(1).eq.channel_rank).and.(id_fromblock(2).eq.3))then
       call saferead_dble(col(3),get_HBresults_obsratio)
      endif
     endif
    enddo
  endif

  deallocate(col)

 end function get_HBresults_obsratio
   
!********************************************************
 function get_HBresults_channel_id(channel_rank)
 !********************************************************
  implicit none
  integer :: id_fromblock(2)
  integer :: blocklines(2)
  integer :: i
  integer, intent(in) :: channel_rank
  integer :: stringlength,stringlength2
  character(len=200) :: temp_line
  integer :: get_HBresults_channel_id
  character(len=50),allocatable :: col(:)

  allocate(col(3))
  get_HBresults_channel_id=0

  blocklines=line_numbers_of_block('HiggsBoundsResults')
  if(minval(blocklines).le.0)then
   write(*,*)'Problem reading the block HiggsBoundsResults'
  else

    do i=1+blocklines(1),blocklines(2)
     temp_line=edited_contents_of_input_file(i)%line

     if(trim(temp_line).ne.'')then

      stringlength =len(temp_line)
      stringlength2=len(col) 

      call split_into_col(temp_line,stringlength,col,stringlength2)

      call saferead_int(col(1),id_fromblock(1))
      call saferead_int(col(2),id_fromblock(2))
      if((id_fromblock(1).eq.channel_rank).and.(id_fromblock(2).eq.1))then
       call saferead_int(col(3),get_HBresults_channel_id)
      endif
     endif
    enddo
  endif

  deallocate(col)

 end function get_HBresults_channel_id
   
 !********************************************************
 subroutine check_validity(is_valid_point)
 !********************************************************
  implicit none
  integer :: error_id,id_fromblock
  integer :: blocklines(2)
  integer :: i,n
  integer :: stringlength,stringlength2
  character(len=200) :: temp_line
  character(len=6) :: blockname
  logical :: is_valid_point
  character(len=50),allocatable :: col(:)

  allocate(col(2))
  is_valid_point=.True.
  error_id=4

  do n=1,2
   select case(n)
   case(1)
    blockname='spinfo'
   case(2)
    blockname='dcinfo'
   case default
    stop'error in function is_valid_point'
   end select
   blocklines=line_numbers_of_block(blockname)
   if(minval(blocklines).le.0)then
    !write(*,*)'Problem reading the block '//trim(adjustl(blockname))
   else

    do i=1+blocklines(1),blocklines(2)
     temp_line=edited_contents_of_input_file(i)%line

     if(trim(temp_line).ne.'')then

      stringlength =len(temp_line)
      stringlength2=len(col) 

      call split_into_col(temp_line,stringlength,col,stringlength2)

      call saferead_int(col(1),id_fromblock)

      if(error_id.eq.id_fromblock)then
       is_valid_point=.False.
      endif
     endif
    enddo
   endif
  enddo

  deallocate(col)

 end subroutine check_validity 
!********************************************************
 function same_particles(array1,array2)
 !********************************************************
  implicit none
  integer,intent(in) :: array1(:)
  integer,intent(in) :: array2(:)

  integer :: record1,record2
  integer :: ub

  integer :: i,j
  logical :: same_particles

  ub=ubound(array1,dim=1)

  if(ub.ne.ubound(array2,dim=1))then
   stop'problem in function same_particles'
  endif

  same_particles=.True.
  do i=1,ub
   record1=0
   record2=0
   do j=1,ub
    if(array1(i).eq.array1(j))then
     record1=record1+1
    endif
    if(array1(i).eq.array2(j))then
     record2=record2+1
    endif
   enddo
   if(record1.ne.record2)then
    same_particles=.False.
    exit 
   endif
  enddo

 end function same_particles

 !********************************************************
 subroutine reorderPDG(array_in,array_out)
 !********************************************************
  implicit none
  integer,intent(in) :: array_in(:)
  integer,allocatable :: array_in_temp(:)
  integer,intent(out) :: array_out(:)
  integer,parameter :: veryhighnumber=1000000000 !should be higher than all the PDG numbers
  integer :: i, minl

  allocate(array_in_temp(lbound(array_in,dim=1):ubound(array_in,dim=1)))

  array_in_temp=array_in
  array_out=0

  do i=lbound(array_in_temp,dim=1),ubound(array_in_temp,dim=1)
    minl=minloc(array_in_temp,dim=1)
    array_out(i)=array_in_temp(minl)
    array_in_temp(minl)=veryhighnumber
  enddo
   
  deallocate(array_in_temp)

 end subroutine reorderPDG
 !********************************************************
 function line_numbers_of_block(name_of_block)
 !********************************************************
  character(len=*),intent(in) ::name_of_block
  integer :: line_numbers_of_block(2)

  line_numbers_of_block=line_numbers_of_block_or_decay('block',name_of_block)
 end function line_numbers_of_block

 !********************************************************
 function line_numbers_of_decay(particlePDGcode)
 !********************************************************
  character(len=7) ::particlePDGcode_string
  integer :: line_numbers_of_decay(2)
  integer,intent(in) :: particlePDGcode

  write(particlePDGcode_string,'(I7)')particlePDGcode

  line_numbers_of_decay=line_numbers_of_block_or_decay('decay',particlePDGcode_string)

 end function line_numbers_of_decay

 !********************************************************
 function line_numbers_of_block_or_decay(typ,name_of_block)
 !********************************************************
 !Returns the line_number of first and last line of block name_of_block
  integer :: n,nmax,m
  integer :: line_numbers_of_block_or_decay(2)
  character(len=*),intent(in) ::name_of_block
  character(len=50) ::name_of_block_lowercase
  character(len=5),intent(in) :: typ
  integer :: stringlength
  character(len=50),allocatable :: col(:)

  allocate(col(2))

  select case(typ)
  case('block')
  case('decay')
  case default
   stop'wrong input to function line_numbers_of_block_or_decay'
  end select

  if(len(name_of_block).gt.len(name_of_block_lowercase))stop'error in function line_numbers_of_block_or_decay (a)'

  line_numbers_of_block_or_decay=0

  stringlength=len(name_of_block)
  name_of_block_lowercase=adjustl(strtolcase(name_of_block,stringlength))
  m=0
  nmax=ubound(names_of_blocks_and_decays,dim=1)
  do n=1,nmax
    if(trim(names_of_blocks_and_decays(n)%line).eq.trim(name_of_block_lowercase))then
      m=m+1
      line_numbers_of_block_or_decay(1)=names_of_blocks_and_decays(n)%id
      if(n.eq.nmax)then
       line_numbers_of_block_or_decay(2)=ubound(edited_contents_of_input_file,dim=1)
      else
       line_numbers_of_block_or_decay(2)=names_of_blocks_and_decays(n+1)%id-1
      endif
    endif
  enddo

  if(m.gt.1)then
   write(*,*)'Warning: there are two '//trim(adjustl(name_of_block))//'in this SLHA file'
  endif


 end function line_numbers_of_block_or_decay
                                                                                                           
 !****************************************************
 function getSLHAfilelength(fileid)
 !****************************************************
 ! calculates file length and checks for errors
 !**************************************************** 
  implicit none
  !--------------------------------------input 
  integer fileid   
  !-----------------------------------function 
  integer :: getSLHAfilelength
  !-----------------------------------internal      
  integer :: n,ios
  character(LEN=5) :: filechar
  character(LEN=20) :: sample
  !-------------------------------------------      
  
  write(filechar,'(I5)')fileid      

  !this will count the number of lines in the file, including the last one
  n = 0                        
  do
   read(fileid,'(a)',iostat=ios) sample

   if(ios.lt.0)then
    exit
   elseif(ios.gt.0) then
    stop 'error in SLHA input file'
   endif
   n = n + 1               
  enddo            

  if(n.eq.0)stop'File is empty'
      
  getSLHAfilelength=n
  rewind(fileid)

 end function getSLHAfilelength

 !****************************************************
end module SLHA_manip
!********************************************************
