! This file is part of HiggsBounds
!  -KW
!******************************************************************
module S95tables_type3
!******************************************************************
 use S95tables_type2
 implicit none

!'checktables' checks that the tables being read in from the 
!text files are in the expected input. 
!Now that we have the binary form of the tables, it's not so important
!to speed up the reading in of the text tables.
!Therefore, might as well leave checktables defined
#define checktables

!oldformat is without the columns Mspecial M0
!#define oldformat 

 !table type 3------------------------------
 type table3
  double precision :: xmax1,xmin1,xmax2,xmin2,sep1,sep2,deltax
  double precision :: zmin,zmax,zsep
  integer :: id,nx1,nx2,nz
  integer :: type_of_assoc_table,id_of_assoc_table 
  logical :: needs_M2_gt_2M1
  double precision, allocatable :: dat(:,:,:)
 end type
 !------------------------------------------

 integer,parameter :: ntable3=10  

 integer,parameter :: skip=3
#ifdef SETCLSBTABLEDIR
 character(len=100),parameter :: pathtofolder&
&=SETCLSBTABLEDIR
#else
 character(len=100),parameter :: pathtofolder='/home/Karina/Work/HiggsBounds/clsb_tables/'
#endif
 double precision,parameter :: deltaMhinfilename=10.0D0 ! results split into sections of size deltaMhinfilename. check carefully if deltaMhinfilename is not a whole number
 integer,parameter :: deltaMhinfilenamebits=5  !which is then split into 5 sections

 integer,parameter :: clsb_t3_dat_3rdcol_clsb       =1
 integer,parameter :: clsb_t3_dat_3rdcol_chisq      =2

#ifdef oldformat
 integer,parameter :: size_of_third_dim_of_dat=2
 integer,parameter :: clsb_t3_dat_3rdcol_Mspecial   =0
 integer,parameter :: clsb_t3_dat_3rdcol_M0         =0
#else
 integer,parameter :: size_of_third_dim_of_dat=4
 integer,parameter :: clsb_t3_dat_3rdcol_Mspecial   =3
 integer,parameter :: clsb_t3_dat_3rdcol_M0         =4
#endif

 character(len=300),allocatable :: filename(:) 
 character(len=100) :: foldername

 !----------------------------------
 type(table3),allocatable :: clsb_t3(:) 
 !----------------------------------

 integer,parameter :: file_id_3=10 !same as file_id_common in usefulbits.f90

 contains

 !************************************************************ 
 subroutine initializetables_type3_blank(tablet3)
 !***********************************************************  
 ! still leaves dat unallocated
  integer:: i
  type(table3) :: tablet3(:)

  do i=lbound(tablet3,dim=1),ubound(tablet3,dim=1)
   tablet3(i)%xmax1       = -1.0D0
   tablet3(i)%xmax2       = -1.0D0
   tablet3(i)%xmin1       = -1.0D0
   tablet3(i)%xmin2       = -1.0D0
   tablet3(i)%sep1        = -1.0D0
   tablet3(i)%sep2        = -1.0D0
   tablet3(i)%deltax      = -1.0D0

   tablet3(i)%zmin    = -1.0D0
   tablet3(i)%zmax    = -1.0D0
   tablet3(i)%zsep    = -1.0D0

   tablet3(i)%id          = -1
   tablet3(i)%nx1         = -1
   tablet3(i)%nx2         = -1
   tablet3(i)%nz          = -1

   tablet3(i)%type_of_assoc_table = -1
   tablet3(i)%id_of_assoc_table   = -1
   tablet3(i)%needs_M2_gt_2M1     = .False.
  enddo
 
 end subroutine initializetables_type3_blank

 !************************************************************ 
 subroutine initializetables3(clsb_t3)
 !***********************************************************  
  implicit none

  !--------------------------------------input
  type(table3) :: clsb_t3(:)
  !-----------------------------------internal
  integer :: x,xbeg,xend
  double precision :: small
  !-------------------------------------------
  small=1.0D-4

  xbeg=lbound(clsb_t3,dim=1)
  xend=ubound(clsb_t3,dim=1)
  
  allocate(filename(xbeg:xend))    
  x=xbeg-1

  x=x+1
  filename(x)       = 'h2h1_bbbb_090307_1_prodah_dechbb'
  clsb_t3(x)%id     = 1800
  clsb_t3(x)%type_of_assoc_table = 2
  clsb_t3(x)%id_of_assoc_table   = 180
  clsb_t3(x)%xmax2  = 180.0D0
  clsb_t3(x)%xmin2  = 1.0D0
  clsb_t3(x)%xmax1  = 180.0D0
  clsb_t3(x)%xmin1  = 1.0D0
  clsb_t3(x)%sep1   = 1.0D0 
  clsb_t3(x)%sep2   = 1.0D0 

  x=x+1
  filename(x)       = 'h2h1_bbtt_090307_1_prodah_dechbt'
  clsb_t3(x)%id     = 2300
  clsb_t3(x)%type_of_assoc_table = 2
  clsb_t3(x)%id_of_assoc_table   = 230
  clsb_t3(x)%xmax2  = 180.0D0
  clsb_t3(x)%xmin2  = 1.0D0
  clsb_t3(x)%xmax1  = 180.0D0
  clsb_t3(x)%xmin1  = 1.0D0
  clsb_t3(x)%sep1   = 1.0D0 
  clsb_t3(x)%sep2   = 1.0D0 

  x=x+1
  filename(x)       = 'h2h1_tatatata_090307_1_prodah_dechtata'
  clsb_t3(x)%id     = 1900
  clsb_t3(x)%type_of_assoc_table = 2
  clsb_t3(x)%id_of_assoc_table   = 190
  clsb_t3(x)%xmax2  = 180.0D0
  clsb_t3(x)%xmin2  = 1.0D0
  clsb_t3(x)%xmax1  = 180.0D0
  clsb_t3(x)%xmin1  = 1.0D0
  clsb_t3(x)%sep1   = 1.0D0 
  clsb_t3(x)%sep2   = 1.0D0 

  x=x+1
  filename(x)       = 'h2h1_ttbb_090307_1_prodah_dechtb'
  clsb_t3(x)%id     = 2400
  clsb_t3(x)%type_of_assoc_table = 2
  clsb_t3(x)%id_of_assoc_table   = 240
  clsb_t3(x)%xmax2  = 180.0D0
  clsb_t3(x)%xmin2  = 1.0D0
  clsb_t3(x)%xmax1  = 180.0D0
  clsb_t3(x)%xmin1  = 1.0D0
  clsb_t3(x)%sep1   = 1.0D0 
  clsb_t3(x)%sep2   = 1.0D0 

  x=x+1
  filename(x)       = 'h2h1_h1h1h1_090307_1_prodah_dechaa'
  clsb_t3(x)%id     = 2000
  clsb_t3(x)%type_of_assoc_table = 2
  clsb_t3(x)%id_of_assoc_table   = 200
  clsb_t3(x)%xmax2  = 180.0D0
  clsb_t3(x)%xmin2  = 1.0D0 ! except there's no numbers for this since Mh2 must be > Mh1
  clsb_t3(x)%xmax1  = 90.0D0
  clsb_t3(x)%xmin1  = 1.0D0
  clsb_t3(x)%sep1   = 1.0D0 
  clsb_t3(x)%sep2   = 1.0D0 

  x=x+1
  filename(x)       = 'h2h1_tttttt_090307_1_prodah_dechtttttt'
  clsb_t3(x)%id     = 2100
  clsb_t3(x)%type_of_assoc_table = 2
  clsb_t3(x)%id_of_assoc_table   = 210
  clsb_t3(x)%xmax2  = 180.0D0
  clsb_t3(x)%xmin2  = 1.0D0 ! except there's no numbers for this since Mh2 must be > Mh1
  clsb_t3(x)%xmax1  = 90.0D0
  clsb_t3(x)%xmin1  = 1.0D0
  clsb_t3(x)%sep1   = 1.0D0 
  clsb_t3(x)%sep2   = 1.0D0
 
  !taken out temporily, because S95_t2 id 220 has also been taken out (see
  !notes for table S95_t2 id 220
  !x=x+1
  !filename(x)       = 'h2z_bbttz_090307_1_prodzh_dechbbtt'
  !clsb_t3(x)%id     = 2200
  !clsb_t3(x)%type_of_assoc_table = 2
  !clsb_t3(x)%id_of_assoc_table   = 220
  !clsb_t3(x)%xmax2  = 180.0D0
  !clsb_t3(x)%xmin2  = 1.0D0 ! except there's no numbers for this since Mh2 must be > Mh1
  !clsb_t3(x)%xmax1  = 90.0D0
  !clsb_t3(x)%xmin1  = 1.0D0
  !clsb_t3(x)%sep1   = 1.0D0 
  !clsb_t3(x)%sep2   = 1.0D0 

  x=x+1
  filename(x)       = 'h2z_h1h1z_090307_1_prodzh_dechaa'
  clsb_t3(x)%id     = 1500
  clsb_t3(x)%type_of_assoc_table = 2
  clsb_t3(x)%id_of_assoc_table   = 150
  clsb_t3(x)%xmax2  = 180.0D0
  clsb_t3(x)%xmin2  = 1.0D0 ! except there's no numbers for this since Mh2 must be > Mh1
  clsb_t3(x)%xmax1  = 90.0D0
  clsb_t3(x)%xmin1  = 1.0D0
  clsb_t3(x)%sep1   = 1.0D0 
  clsb_t3(x)%sep2   = 1.0D0 

  x=x+1
  filename(x)       = 'h2z_ttttz_090307_1_prodzh_dechtttt'
  clsb_t3(x)%id     = 1600
  clsb_t3(x)%type_of_assoc_table = 2
  clsb_t3(x)%id_of_assoc_table   = 160
  clsb_t3(x)%xmax2  = 180.0D0
  clsb_t3(x)%xmin2  = 1.0D0 ! except there's no numbers for this since Mh2 must be > Mh1
  clsb_t3(x)%xmax1  = 90.0D0
  clsb_t3(x)%xmin1  = 1.0D0
  clsb_t3(x)%sep1   = 1.0D0 
  clsb_t3(x)%sep2   = 1.0D0 

  x=x+1
  filename(x)       = 'h2z_bbz_090307_1_prodzh_dechbb'
  clsb_t3(x)%id     = 1420
  clsb_t3(x)%type_of_assoc_table = 1
  clsb_t3(x)%id_of_assoc_table   = 142
  clsb_t3(x)%xmax1  = 180.0D0
  clsb_t3(x)%xmin1  = 0.1D0
  clsb_t3(x)%xmax2  = 20.0D0
  clsb_t3(x)%xmin2  = 20.0D0
  clsb_t3(x)%sep2   = 1.0D0 !doesn't matter what this value is as long as it is not zero
  clsb_t3(x)%sep1   = 0.1D0 

  x=x+1
  filename(x)       = 'h2z_tataz_090307_1_prodzh_dechtata'
  clsb_t3(x)%id     = 1430
  clsb_t3(x)%type_of_assoc_table = 1
  clsb_t3(x)%id_of_assoc_table   = 143
  clsb_t3(x)%xmax1  = 180.0D0
  clsb_t3(x)%xmin1  = 0.1D0
  clsb_t3(x)%xmax2  = 20.0D0
  clsb_t3(x)%xmin2  = 20.0D0
  clsb_t3(x)%sep2   = 1.0D0 !doesn't matter what this value is as long as it is not zero
  clsb_t3(x)%sep1   = 0.1D0 

  ! checks we've filled the whole array
  if(x.ne.xend)then
   stop'error in initializetables3 (a)'
  endif
 
  do x=xbeg,xend
   clsb_t3(x)%zmin    = -3.0D0
   clsb_t3(x)%zmax    = 0.03D0
   clsb_t3(x)%zsep    = 0.03D0

   clsb_t3(x)%nx2     = nint((clsb_t3(x)%xmax2-clsb_t3(x)%xmin2)/clsb_t3(x)%sep2)+1
   clsb_t3(x)%nx1     = nint((clsb_t3(x)%xmax1-clsb_t3(x)%xmin1)/clsb_t3(x)%sep1)+1
   clsb_t3(x)%nz      = nint((clsb_t3(x)%zmax -clsb_t3(x)%zmin )/clsb_t3(x)%zsep)+1

   if(abs(clsb_t3(x)%sep1).lt.small)stop'problem in subroutine initializetables3 b1'
   if(abs(clsb_t3(x)%sep2).lt.small)stop'problem in subroutine initializetables3 b2'
   if(abs(clsb_t3(x)%zsep).lt.small)stop'problem in subroutine initializetables3 b3'
   if(clsb_t3(x)%xmax1.lt.clsb_t3(x)%xmin1)stop'problem in subroutine initializetables3 b4'
   if(clsb_t3(x)%xmax2.lt.clsb_t3(x)%xmin2)stop'problem in subroutine initializetables3 b5'
   if(clsb_t3(x)%zmax .lt.clsb_t3(x)%zmin )stop'problem in subroutine initializetables3 b6'

   !problem: can't set clsb_t3(x)%needs_M2_gt_2M1 here because it needs 
   !S95_t2
   !if(clsb_t3(x)%type_of_assoc_table.eq.2)then
   !  table2num=t2elementnumberfromid(S95_t2,clsb_t3(x)%id_of_assoc_table)
   !  clsb_t3(x)%needs_M2_gt_2M1 =  S95_t2(table2num)%needs_M2_gt_2M1
   !endif
   !instead, call fillt3needs_M2_gt_2M1

  enddo

 end subroutine initializetables3
 !*********************************************************** 
 subroutine fillt3needs_M2_gt_2M1(t3,t2)
  !--------------------------------------input
  type(table3) :: t3(:)
  type(table2), intent(in) :: t2(:)
  !-----------------------------------internal
  integer :: n,x
  !-------------------------------------------

  do x=lbound(t3,dim=1),ubound(t3,dim=1)
   if(t3(x)%type_of_assoc_table.eq.2)then
     n=t2elementnumberfromid(t2,t3(x)%id_of_assoc_table)
     t3(x)%needs_M2_gt_2M1 = t2(n)%needs_M2_gt_2M1
   endif
  enddo

 end subroutine fillt3needs_M2_gt_2M1
 !*********************************************************** 
 function t3elementnumberfromid(t3,id)
  !--------------------------------------input
  type(table3), intent(in) :: t3(:)
  integer, intent(in) :: id
  !-----------------------------------function
  integer :: t3elementnumberfromid
  !-----------------------------------internal
  integer :: n,x
  !-------------------------------------------

  n=0
  do x=lbound(t3,dim=1),ubound(t3,dim=1)
   if(t3(x)%id.eq.id)then
    n=n+1
    t3elementnumberfromid=x
   endif
  enddo

  if(n.ne.1)stop'problem in function t3elementnumberfromid 1'

 end function t3elementnumberfromid

 !*********************************************************** 
 subroutine fill_t1dat_from_t3(t1,t3,zi,ftype)
 ! note that t1%sep,t1%xmin,t1%xmax,t1%nx should already be set
 ! and t1%dat should be allocated
 !*********************************************************** 
  !--------------------------------------input
  type(table3), intent(in) :: t3
  integer, intent(in) :: zi,ftype
  !-----------------------------------internal
  integer :: a,ji,jilow,jihigh
  !integer :: j,i,jtemp,itemp,itot
  !------------------------------------changed
  type(table1) :: t1
  !-------------------------------------------

  if(    t1%sep.eq.0)then
    stop'problem in fill_t1dat_from_t3 (2a)'
  elseif(t1%nx.eq.0)then
    stop'problem in fill_t1dat_from_t3 (2b)'
  elseif(nint((t1%xmax-t1%xmin)/t1%sep)+1.ne.t1%nx)then
    stop'problem in fill_t1dat_from_t3 (2c)'
  endif

  if(t1%sep.ne.t3%sep1)stop'problem in fill_t1dat_from_t3 (2d)'

  if(ftype.lt.lbound(t3%dat,dim=3))stop'problem in fill_t1dat_from_t3 (2a)'
  if(ftype.gt.ubound(t3%dat,dim=3))stop'problem in fill_t1dat_from_t3 (2b)'

  if(zi.lt.lbound(t3%dat,dim=2))stop'problem in fill_t1dat_from_t3 (3a)'
  if(zi.gt.ubound(t3%dat,dim=2))stop'problem in fill_t1dat_from_t3 (3b)'
  
  if(ubound(t1%dat,dim=1).ne.t1%nx)stop'problem in fill_t1dat_from_t3 (5a)'

  select case(t3%type_of_assoc_table)
  case(1)
    if(t1%xmin.lt.t3%xmin1)stop'problem in fill_t1dat_from_t3 (4a)'
    if(t1%xmax.gt.t3%xmax1)stop'problem in fill_t1dat_from_t3 (4b)'


    ! at the moment, want t1%xmin, t1%xmax to correspond to points in t3
    ! can edit this later if needed
    if(abs((t1%xmin-t3%xmin1)/t3%sep1-dble(nint((t1%xmin-t3%xmin1)/t3%sep1))).gt.1.0D-3)then
     stop'problem in fill_t1dat_from_t3 (6a)'
    endif
    if(abs((t1%xmax-t3%xmax1)/t3%sep1-dble(nint((t1%xmax-t3%xmax1)/t3%sep1))).gt.1.0D-3)then
     stop'problem in fill_t1dat_from_t3 (6b)'
    endif
  case(2)
    if(t1%xmin.lt.t3%xmin2)stop'problem in fill_t1dat_from_t3 (5a)'
    if(t1%xmax.gt.t3%xmax2)stop'problem in fill_t1dat_from_t3 (5b)'

    if(abs((t1%xmin-t3%xmin2)/t3%sep2-dble(nint((t1%xmin-t3%xmin2)/t3%sep2))).gt.1.0D-3)then
     stop'problem in fill_t1dat_from_t3 (5c)'
    endif
    if(abs((t1%xmax-t3%xmax2)/t3%sep2-dble(nint((t1%xmax-t3%xmax2)/t3%sep2))).gt.1.0D-3)then
     stop'problem in fill_t1dat_from_t3 (5d)'
    endif
  case default
    stop'problem in fill_t1dat_from_t3 (1a)'
  end select 


  t1%dat = -1.0D0

  select case(t3%type_of_assoc_table)
  case(1)

   a=0

   !t3%nx2=1 for t3%type_of_assoc_table=1
   jilow =nint((t1%xmin-t3%xmin1)/t3%sep1)+1
   jihigh=nint((t1%xmax-t3%xmin1)/t3%sep1)+1

   do ji=jilow,jihigh
    a=a+1
    t1%dat(a,1)=t3%dat(ji,zi,ftype)
   enddo

   if(a.ne.ubound(t1%dat,dim=1))stop'problem in fill_t1dat_from_t3 (10)'

  case(2)
   !This would need extra arguments to the function, since one
   !of the masses would need to be kept constant
   stop'problem in subroutine fill_t1dat_from_t3'
  case default
    stop'problem in fill_t1dat_from_t3 (1)'
  end select 

 end subroutine fill_t1dat_from_t3
 !*********************************************************** 
 subroutine fill_slices_t2_from_slices_of_t3(t3,zi_selection,ftype_selection,slices_t2)
 ! if this subroutine is used,  
 ! don't forget to deallocate slices_t2(x)%dat at some point
 !*********************************************************** 
 implicit none 
  !--------------------------------------input
  type(table3), intent(in) :: t3
  type(table2), intent(out) :: slices_t2(:)
  integer, intent(in) :: zi_selection(:)
  integer, intent(in) :: ftype_selection(:)
  !-----------------------------------internal
  integer :: i,j,k,x,ji
  integer :: itot
  !-------------------------------------------

  if(lbound(zi_selection,dim=1).ne.lbound(slices_t2,dim=1))then
    stop'problem in fill_slices_t2_from_slices_of_t3 1a'
  endif
  if(ubound(zi_selection,dim=1).ne.ubound(slices_t2,dim=1))then
    stop'problem in fill_slices_t2_from_slices_of_t3 1b'
  endif

  do x=lbound(ftype_selection,dim=1),ubound(ftype_selection,dim=1)
   if(ftype_selection(x).lt.lbound(t3%dat,dim=3))stop'problem in fill_slices_t2_from_slices_of_t3 3a'
   if(ftype_selection(x).gt.ubound(t3%dat,dim=3))stop'problem in fill_slices_t2_from_slices_of_t3 3b'
  enddo

  do x=lbound(slices_t2,dim=1),ubound(slices_t2,dim=1)

   if(zi_selection(x).lt.lbound(t3%dat,dim=2))stop'problem in fill_slices_t2_from_slices_of_t3 4a'
   if(zi_selection(x).gt.ubound(t3%dat,dim=2))stop'problem in fill_slices_t2_from_slices_of_t3 4b'


   slices_t2(x)%id          =  t3%id    
   slices_t2(x)%nx1         =  t3%nx1        
   slices_t2(x)%nx2         =  t3%nx2  
   slices_t2(x)%xmax1       =  t3%xmax1  
   slices_t2(x)%xmax2       =  t3%xmax2
   slices_t2(x)%xmin1       =  t3%xmin1
   slices_t2(x)%xmin2       =  t3%xmin2
   slices_t2(x)%sep1        =  t3%sep1 
   slices_t2(x)%sep2        =  t3%sep2
   slices_t2(x)%deltax      =  t3%deltax

   slices_t2(x)%z      =  t3%zmin + dble(zi_selection(x)-1)*t3%zsep

   slices_t2(x)%needs_M2_gt_2M1      =  t3%needs_M2_gt_2M1

   allocate( slices_t2(x)%dat(slices_t2(x)%nx2,slices_t2(x)%nx1, &
           & lbound(ftype_selection,dim=1):ubound(ftype_selection,dim=1))) 
   slices_t2(x)%dat = -1.0D0

   ji=0
   do j=1,slices_t2(x)%nx2

    select case(t3%type_of_assoc_table)
    case(2)!how general is this?
     itot=int(dble(j)*dble(t3%nx1)/dble(t3%nx2))
    case(0,1)
     itot=t3%nx1
    case default
     stop'problem in fill_slices_t2_from_slices_of_t3 5'
    end select

    do i=1,itot
     ji=ji+1
     do k=lbound(ftype_selection,dim=1),ubound(ftype_selection,dim=1)
       slices_t2(x)%dat(j,i,k)=t3%dat(ji,zi_selection(x),ftype_selection(k))
     enddo
    enddo
   enddo 

   if(ji.ne.size(t3%dat,dim=1))stop'problem in fill_slices_t2_from_slices_of_t3 6'

  enddo

 end subroutine fill_slices_t2_from_slices_of_t3
 !*********************************************************** 
 subroutine readclsbfiles_binary
 !*********************************************************** 
  implicit none
  integer :: z
  integer :: n_xy_combinations

  do z=lbound(clsb_t3,dim=1),ubound(clsb_t3,dim=1)

   n_xy_combinations=clsb_t3(z)%nx2*(clsb_t3(z)%nx1+int(clsb_t3(z)%nx1/clsb_t3(z)%nx2))/2

   allocate(clsb_t3(z)%dat(n_xy_combinations,clsb_t3(z)%nz,size_of_third_dim_of_dat))
  
   foldername='csboutput_trans_binary/csboutput_outfile_fullclsb_ee_'//trim(adjustl(filename(z)))
   open(file_id_3,file = trim(adjustl(pathtofolder)) // &
      &  trim(adjustl(foldername))//'.binary',form='unformatted')

   read(file_id_3)clsb_t3(z)%dat

   close(file_id_3)
   write(*,*)'finished reading in ',trim(filename(z))
  enddo

 end subroutine readclsbfiles_binary
 !*********************************************************** 
 subroutine writeclsbfiles_binary
 !*********************************************************** 
 implicit none
 integer :: z

  do z=lbound(clsb_t3,dim=1),ubound(clsb_t3,dim=1)

   foldername='csboutput_trans_binary/csboutput_outfile_fullclsb_ee_'//trim(adjustl(filename(z)))
   open(file_id_3,file = trim(adjustl(pathtofolder)) // &
      &  trim(adjustl(foldername))//'.binary',form='unformatted')

   write(file_id_3)clsb_t3(z)%dat

   close(file_id_3)
  enddo

 end subroutine writeclsbfiles_binary
 !*********************************************************** 
 subroutine readclsbfiles(clsb_t3_id)
 !*********************************************************** 
 use usefulbits, only : vsmall
 implicit none
  !--------------------------------------input
  integer,intent(in) :: clsb_t3_id
  !-----------------------------------internal
  double precision :: Mh2_in, Mh1_in, log10sf_in, clsb_in, chisq_in,Mspecial_in,M0_in
  double precision :: Mhinfilename_initial,Mhinfilename_xmin
  integer :: i,a,n_xy_combinations,b,c
  integer :: n,j,x,n_csbfile,y,jmax,ymax
  integer, allocatable :: imax(:)
  integer :: Mhbeg,Mhend
  character(len=7)   :: Mhstringpart1
  character(len=2)   :: Mhstringpart2
  character(len=3)   :: Mhbegstring,Mhendstring
  character(len=1)   :: Mhdigit
  integer ::  z
  !------------------------------------------

   n=0
   do x=lbound(clsb_t3,dim=1),ubound(clsb_t3,dim=1)
    if(clsb_t3(x)%id.eq.clsb_t3_id)then
     n=n+1
     z=x
    endif
   enddo

   if(n.ne.1)then
    write(*,*)'n=',n
    stop'problem in readclsbfiles 1'
   endif

   a=0 !a labels the number of data points in entire table (spread over many files)
   b=0 !b labels the number of (x,y) points i.e. will be compared to n_xy_combinations

   n_xy_combinations=clsb_t3(z)%nx2*(clsb_t3(z)%nx1+int(clsb_t3(z)%nx1/clsb_t3(z)%nx2))/2

   allocate(clsb_t3(z)%dat(n_xy_combinations,clsb_t3(z)%nz,size_of_third_dim_of_dat))

#ifdef oldformat
   foldername='csboutput_trans/csboutput_outfile_fullclsb_ee_'//trim(adjustl(filename(z)))
#else
   foldername='csboutput_trans_withMexclM0/csboutput_outfile_fullclsb_ee_'//trim(adjustl(filename(z)))
#endif
  
   ! n.b. this is not (clsb_t3(z)%xmax2-clsb_t3(z)%xmin2)/deltaMhinfilename
   ! because of the file format: the first file is shorter than would be the case if the above expression 
   ! was correct

   if(abs(deltaMhinfilename-nint(deltaMhinfilename)).gt.vsmall)stop'deltaMhinfilename must be approx an integer'

   select case(clsb_t3(z)%type_of_assoc_table)
   case(1)
    ymax=nint((clsb_t3(z)%xmax1)/deltaMhinfilename) 
    Mhinfilename_xmin=clsb_t3(z)%xmin1
    jmax=1
   case(2)
    ymax=nint((clsb_t3(z)%xmax2)/deltaMhinfilename) 
    Mhinfilename_xmin=clsb_t3(z)%xmin2
    jmax=nint(deltaMhinfilename/dble(deltaMhinfilenamebits)/clsb_t3(z)%sep2)
   case default
    stop'error in subroutine readclsbfiles'
   end select

   allocate(imax(jmax))

   do y=1,ymax
    do x=1,deltaMhinfilenamebits    
     
     Mhinfilename_initial=Mhinfilename_xmin + dble(y-1)*deltaMhinfilename + dble(x-1)*deltaMhinfilename/dble(deltaMhinfilenamebits)

     select case(clsb_t3(z)%type_of_assoc_table)
     case(1)
      imax=nint(deltaMhinfilename/dble(deltaMhinfilenamebits)/clsb_t3(z)%sep1)
     case(2)
      do j=1,jmax
       imax(j)=int((Mhinfilename_initial+clsb_t3(z)%sep2*(j-1))*(clsb_t3(z)%xmax1/clsb_t3(z)%xmax2)/clsb_t3(z)%sep2)
      enddo
     case default
      stop'error in subroutine readclsbfiles'
     end select
  
     Mhbeg=  ( int( (Mhinfilename_initial)/deltaMhinfilename ) )*nint(deltaMhinfilename) 
     Mhend= Mhbeg+10
    
     write(Mhbegstring,'(I3)')Mhbeg
     write(Mhendstring,'(I3)')Mhend

     select case(Mhbeg)
     case(0:9)
       Mhbegstring='00'//adjustl(Mhbegstring)
     case(10:99)
       Mhbegstring='0'//adjustl(Mhbegstring)
     case default
     end select
  
     select case(Mhend)
     case(0:9)
       Mhendstring='00'//adjustl(Mhendstring)
     case(10:99)
       Mhendstring='0'//adjustl(Mhendstring)
     case default
     end select

     Mhstringpart1=Mhbegstring//'_'//Mhendstring
     !write(*,*)'~'//Mhstringpart1//'~'

     if(x.gt.9)stop'have not done this case yet'
     write(Mhdigit,'(I1)')x
     Mhstringpart2='_'//Mhdigit

     open(file_id_3,file = trim(adjustl(pathtofolder)) // &
      &  trim(adjustl(foldername))//'_loopmh'//Mhstringpart1//Mhstringpart2//'.txt')

     c=0 ! c labels the number of lines read in from a single file
  
     do n=1,skip
      read(file_id_3,*)
     enddo

     n_csbfile=skip
     do j=1,jmax

#ifdef checktables
      do i= 1, imax(j)
       do n=1,clsb_t3(z)%nz
        n_csbfile=n_csbfile+1
       enddo
      enddo
#endif

      do i=1,imax(j)
       b=b+1
       do n=1,clsb_t3(z)%nz
        c=c+1
        a=a+1

#ifdef oldformat
        read(file_id_3,*)Mh2_in, Mh1_in, log10sf_in, clsb_in, chisq_in
#else
        read(file_id_3,*)Mh2_in, Mh1_in, log10sf_in, clsb_in, chisq_in, Mspecial_in, M0_in
#endif

#ifdef checktables
        !check Mh1_in,Mh2_in are as expected
        select case(clsb_t3(z)%type_of_assoc_table)
        case(1)
          if(abs(Mh1_in-( Mhinfilename_initial + dble(i-1)*clsb_t3(z)%sep1)) .gt.0.0001D0 )stop'error: Mh1_in'
          if(abs(Mh2_in-(clsb_t3(z)%xmin2 + dble(j-1)*clsb_t3(z)%sep2)) .gt.0.0001D0 )     stop'error: Mh2_in'
        case(2)
          if(abs(Mh1_in-(clsb_t3(z)%xmin1 + dble(i-1)*clsb_t3(z)%sep1)) .gt.vsmall   )     stop'error: Mh1_in'
          if(abs(Mh2_in-( Mhinfilename_initial + dble(j-1)*clsb_t3(z)%sep2)) .gt.0.0001D0 )stop'error: Mh2_in'
        case default
          stop'error in subroutine readclsbfiles'
        end select

        !check log10sf_in is as expected
        if(  abs(log10sf_in-(clsb_t3(z)%zmin+clsb_t3(z)%zsep*dble(n-1))) .gt.0.0001D0 )stop'error: log10sf_in'
#endif

        clsb_t3(z)%dat(b,n,clsb_t3_dat_3rdcol_clsb    )=clsb_in
        clsb_t3(z)%dat(b,n,clsb_t3_dat_3rdcol_chisq   )=chisq_in 


#ifndef oldformat
        clsb_t3(z)%dat(b,n,clsb_t3_dat_3rdcol_Mspecial    )=Mspecial_in
        clsb_t3(z)%dat(b,n,clsb_t3_dat_3rdcol_M0          )=M0_in 
#endif

       enddo
      enddo
     enddo 

     close(file_id_3)

#ifdef checktables
     if((n_csbfile-skip).ne.c)then
       write(*,*)'n_csbfile,skip=',n_csbfile,skip
       write(*,*)'c=',c
       stop'error: n_csbfile,skip'
     endif 
#endif
    enddo
   enddo

   deallocate(imax)

#ifdef checktables
   if(abs(clsb_t3(z)%xmax2-Mh2_in).gt.vsmall)stop'error in clsb_t3(z)%xmax2'
   if(abs(clsb_t3(z)%xmax1-Mh1_in).gt.vsmall)then
    write(*,*)'hello clsb_t3(z)%xmax1,Mh1_in',clsb_t3(z)%xmax1,Mh1_in
    stop'error in clsb_t3(z)%xmax1'
   endif
   if(abs(clsb_t3(z)%zmax-log10sf_in).gt.0.001)stop'error in clsb_t3(z)%zmax'
   if((n_xy_combinations).ne.b)stop'error in n_xy_combinations'
   if((n_xy_combinations*clsb_t3(z)%nz).ne.a)stop'error in n_xy_combinations*clsb_t3(z)%nz'
#endif

   write(*,*)'finished reading in ',trim(filename(z))
   call flush(6)

 end subroutine readclsbfiles
 !*********************************************************** 
 function clsb_t3elementnumber_from_S95table(ttype,id)
  implicit none 
  integer,intent(in) :: ttype,id
  integer :: z,n
  integer :: clsb_t3elementnumber_from_S95table

  if(.not.(allocated(clsb_t3)))stop'error in function clsb_t3elementnumber_from_S95table'

  clsb_t3elementnumber_from_S95table=-1

  n=0
  do z=lbound(clsb_t3,dim=1),ubound(clsb_t3,dim=1)
   if((clsb_t3(z)%type_of_assoc_table.eq.ttype).and. &
    & (clsb_t3(z)%id_of_assoc_table  .eq.id   ))then
     n=n+1
     clsb_t3elementnumber_from_S95table=z
   endif
  enddo

  if(n.gt.1)stop'error in function clsb_t3elementnumber_from_S95table'
 
 end function clsb_t3elementnumber_from_S95table
 !*********************************************************** 
 function gett3dat(intable,t3, Mi,Mj,log10sf,selec)
  implicit none
  type(table3), intent(in) :: t3
  double precision, intent(in) :: Mj,Mi,log10sf
  integer :: i,j,z,ji,jtemp,itemp,selec,itot
  logical, intent(out):: intable
  double precision :: ibit,jbit,zbit
  double precision :: gett3dat

  intable=.True.

  i=nint(  (Mi-t3%xmin1)   /t3%sep1)+1
  j=nint(  (Mj-t3%xmin2)   /t3%sep2)+1
  z=nint((log10sf-t3%zmin )/t3%zsep)+1

  if(min(i,j,z).lt.1)then
    intable=.False.
  elseif(i.gt.t3%nx1)then
    intable=.False.
  elseif(j.gt.t3%nx2)then
    intable=.False.
  elseif(z.gt.t3%nz )then
    intable=.False.
  endif

  ibit=(Mi     -(dble(i-1)*t3%sep1+t3%xmin1))/t3%sep1
  jbit=(Mj     -(dble(j-1)*t3%sep2+t3%xmin2))/t3%sep2
  zbit=(log10sf-(dble(z-1)*t3%zsep+t3%zmin ))/t3%zsep

  if(max(ibit,jbit,zbit).gt.0.01D0)then ! Mj,Mi,log10sf are not on a data point
    stop'bad input to function gett3dat'
  endif

  ji=0
  do jtemp=1,j
   if(jtemp.eq.j)then
    itot=i
   else
    select case(t3%type_of_assoc_table)
    case(2)!how general is this?
     itot=int(dble(jtemp)*dble(t3%nx1)/dble(t3%nx2))
    case(0,1)
     itot=t3%nx1
    case default
     stop'problem in function gett3dat'
    end select
   endif
 
   do itemp=1,itot
    ji=ji+1
   enddo
  enddo 
  
  if(intable)then
   gett3dat=t3%dat(ji,z,selec)
  else
   gett3dat=-1.0D0
  endif

 end function gett3dat

end module S95tables_type3
!****************************************************************** _
