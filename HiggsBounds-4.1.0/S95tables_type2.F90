! This file is part of HiggsBounds
!  -KW
!******************************************************************
module S95tables_type2
!******************************************************************
 use S95tables_type1
 implicit none
!#define fast
 
 !table type 2------------------------------
 type table2
  logical :: needs_M2_gt_2M1
  integer :: id,nx1,nx2,particle_x1,particle_x2 !see usefulbits.f90 for key to particle codes n.b. they're NOT pdg
  character(LEN=45) :: label
  character(LEN=3) :: expt      
  double precision :: lumi, energy  
  double precision :: xmax1,xmin1,xmax2,xmin2,sep1,sep2,deltax
  double precision, allocatable :: dat(:,:,:) !in dat(a,b,1:2)...obs=1,pred=2. 1st component of dat is y, 2nd is x 
  double precision :: maxdatval    
  double precision :: z !only used in slices_t2
 end type

 integer,parameter :: file_id_2_exp=10  !same as file_id_common in usefulbits.f90
 integer,parameter :: file_id_2_obs=11
                                    
 !------------------------------------------ 

 contains


 !************************************************************ 
 subroutine initializetables_type2_blank(tablet2)
 !***********************************************************  
 ! still leaves dat unallocated
  integer:: i
  type(table2) :: tablet2(:)

  do i=lbound(tablet2,dim=1),ubound(tablet2,dim=1)
   tablet2(i)%id          = -1
   tablet2(i)%nx1         = -1
   tablet2(i)%nx2         = -1
   tablet2(i)%particle_x1 = -1
   tablet2(i)%particle_x2 = -1
   tablet2(i)%label       = ''
   tablet2(i)%expt        = ''
   tablet2(i)%lumi       = -1.0D0
   tablet2(i)%energy     = -1.0D0   
   tablet2(i)%xmax1       = -1.0D0
   tablet2(i)%xmax2       = -1.0D0
   tablet2(i)%xmin1       = -1.0D0
   tablet2(i)%xmin2       = -1.0D0
   tablet2(i)%sep1        = -1.0D0
   tablet2(i)%sep2        = -1.0D0
   tablet2(i)%deltax      = -1.0D0
   tablet2(i)%maxdatval   = -1.0D0

   tablet2(i)%z   = -1.0D9 !only used in slices_t2

   tablet2(i)%needs_M2_gt_2M1     = .False.
  enddo
 
 end subroutine initializetables_type2_blank

 !*********************************************************** 
 subroutine initializetables2(S95_t2)
 !***********************************************************  
 ! fills S95_t2
 !***********************************************************  
  use store_pathname      
  use usefulbits, only: Hneut,Chineut,Chiplus,small,file_id_common2,not_a_particle
  implicit none

  !--------------------------------------input
  type(table2) :: S95_t2(:)
  !-----------------------------------internal
  integer :: i,tno,j,x,xbeg,xend,k,ios
  character(LEN=2) :: tableno            
  character(len=100),allocatable :: filename(:) 
  double precision :: dummy
  double precision, allocatable :: testrow(:)
  integer :: file_id_arr(2)
  double precision :: maxdatval
  !-------------------------------------------       
  file_id_arr(1)=file_id_2_exp
  file_id_arr(2)=file_id_2_obs

  xbeg=lbound(S95_t2,dim=1)
  xend=ubound(S95_t2,dim=1)  
 
  allocate(filename(xbeg:xend))
  x=xbeg-1
   
  tno=14
  do i=1,8
   x=x+1
   tno=tno+1
   if((x.eq.3).or.(x.eq.7))tno=tno+1 
   write(tableno,'(I2)')tno      
   
   S95_t2(x)%id=tno*10
   S95_t2(x)%expt='LEP'
   S95_t2(x)%energy=0.208D0        
   S95_t2(x)%deltax=0.0D0 
   S95_t2(x)%particle_x1=Hneut
   S95_t2(x)%particle_x2=Hneut

   select case(S95_t2(x)%id)
   case(220,230,240)
    S95_t2(x)%label='hep-ex/0602042 (LEP)'
   case default
    S95_t2(x)%label='hep-ex/0602042, table '//tableno//' (LEP)'
   end select

   S95_t2(x)%sep1=1.0D0
   S95_t2(x)%sep2=1.0D0
   S95_t2(x)%maxdatval   = 1.0D2   
   !S95_t2(x)%OBid=x+2
       
   select case(S95_t2(x)%id)
   case(150,160,220)
    S95_t2(x)%xmin1=1.0D0
    S95_t2(x)%xmax1=60.0D0
    S95_t2(x)%xmin2=2.0D0
    S95_t2(x)%xmax2=120.0D0
    S95_t2(x)%needs_M2_gt_2M1=.True.
   case(180,190,230,240)   
    S95_t2(x)%xmin1=1.0D0
    S95_t2(x)%xmax1=180.0D0
    S95_t2(x)%xmin2=1.0D0
    S95_t2(x)%xmax2=180.0D0 
    S95_t2(x)%needs_M2_gt_2M1=.False.
   case(200,210) 
    S95_t2(x)%xmin1=1.0D0
    S95_t2(x)%xmax1=90.0D0
    S95_t2(x)%xmin2=2.0D0
    S95_t2(x)%xmax2=180.0D0 
    S95_t2(x)%needs_M2_gt_2M1=.True.                           
   case default
   write(*,*)'error in initializetables2 (a)' 
    stop
   end select
   
   filename(x)='table'//tableno//'full' 
  enddo
   
  do i=5,10
   x=x+1
   tno=i
   write(tableno,'(I2)')tno      
   
   S95_t2(x)%id=900+tno
   S95_t2(x)%expt='LEP' 
   S95_t2(x)%energy=0.208D0     
   S95_t2(x)%deltax=0.0D0    
   S95_t2(x)%label='hep-ex/0401026, fig '//trim(adjustl(tableno))//' (OPAL)'
   S95_t2(x)%sep1=1.0D0
   S95_t2(x)%sep2=1.0D0
   S95_t2(x)%maxdatval   = 1.0D6 !these tables are in fb   
       
   select case(tno)
   case(5,6,7,8)
    S95_t2(x)%xmin1=0.0D0
    S95_t2(x)%xmax1=100.0D0
    S95_t2(x)%xmin2=75.0D0
    S95_t2(x)%xmax2=120.0D0
    S95_t2(x)%needs_M2_gt_2M1=.False. 
    S95_t2(x)%particle_x1=Chineut
    S95_t2(x)%particle_x2=Chiplus
   case(9,10)   
    S95_t2(x)%xmin1=0.0D0
    S95_t2(x)%xmax1=100.0D0
    S95_t2(x)%xmin2=50.0D0
    S95_t2(x)%xmax2=200.0D0 
    S95_t2(x)%needs_M2_gt_2M1=.False. 
    S95_t2(x)%particle_x1=Chineut
    S95_t2(x)%particle_x2=Chineut                    
   case default
   write(*,*)'error in initializetables2 (b)' 
    stop
   end select
   
   filename(x)='1026_fig'//trim(adjustl(tableno)) 
  enddo  

  x=x+1    
  S95_t2(x)%id=3381
  S95_t2(x)%expt=' D0'    
  S95_t2(x)%energy=1.96D0  
  S95_t2(x)%deltax=0.0D0   
  S95_t2(x)%label='[hep-ex] arXiv:0905.3381, table I (D0)'
  S95_t2(x)%sep1=0.1D0
  S95_t2(x)%sep2=5.0D0
  S95_t2(x)%maxdatval=1.0D6
  S95_t2(x)%xmin1=0.2D0
  S95_t2(x)%xmax1=3.0D0
  S95_t2(x)%xmin2=80.0D0
  S95_t2(x)%xmax2=200.0D0
  S95_t2(x)%needs_M2_gt_2M1=.True. 
  S95_t2(x)%particle_x1=Hneut
  S95_t2(x)%particle_x2=Hneut
  filename(x)="D0_h-aa-mumumumu_3381" 

  x=x+1    
  S95_t2(x)%id=3382
  S95_t2(x)%expt=' D0'   
  S95_t2(x)%energy=1.96D0    
  S95_t2(x)%deltax=0.0D0      
  S95_t2(x)%label='[hep-ex] arXiv:0905.3381, table II (D0)'
  S95_t2(x)%sep1=0.2D0
  S95_t2(x)%sep2=5.0D0
  S95_t2(x)%maxdatval=1.0D6
  S95_t2(x)%xmin1=3.6D0
  S95_t2(x)%xmax1=19.0D0
  S95_t2(x)%xmin2=85.0D0
  S95_t2(x)%xmax2=200.0D0
  S95_t2(x)%needs_M2_gt_2M1=.True. 
  S95_t2(x)%particle_x1=Hneut
  S95_t2(x)%particle_x2=Hneut
  filename(x)="D0_h-aa-tautaumumu_3381"

  x=x+1 
  S95_t2(x)%id=6227
  S95_t2(x)%expt=' D0'
  S95_t2(x)%label='D0 Note 6227'
  S95_t2(x)%energy=1.96D0    
  S95_t2(x)%sep1=0.04D0 
  S95_t2(x)%sep2=10.0D0 
  S95_t2(x)%maxdatval=1.0D6
  S95_t2(x)%xmin1=0.06D0
  S95_t2(x)%xmax1=0.18D0
  S95_t2(x)%xmin2=90.0D0
  S95_t2(x)%xmax2=300.0D0
  S95_t2(x)%needs_M2_gt_2M1=.False. 
  S95_t2(x)%particle_x1=not_a_particle
  S95_t2(x)%particle_x2=Hneut
  filename(x)='D0_h-bb_h-tautau_comb_5.2-7.3fb_6227'  

  ! checks we've filled the whole array  
  if(x.ne.xend)then
   write(*,*)'error in initializetables2 (c)',x,xend
   stop
  endif  
  
  ! read in the tables
  do x=xbeg,xend
   
   S95_t2(x)%nx2=nint((S95_t2(x)%xmax2-S95_t2(x)%xmin2)/S95_t2(x)%sep2)+1
   S95_t2(x)%nx1=nint((S95_t2(x)%xmax1-S95_t2(x)%xmin1)/S95_t2(x)%sep1)+1             

   allocate(S95_t2(x)%dat(S95_t2(x)%nx2,S95_t2(x)%nx1,2)) 
  enddo

  ! read in the tables
  open(file_id_common2,file = trim(adjustl(pathname))//'Expt_tables/' // &
      &  'S95_t2.binary',form='unformatted')

  read(file_id_common2,iostat=ios)S95_t2(xbeg)%dat
  if(ios.eq.0)then
    do x=xbeg+1,xend
     read(file_id_common2)S95_t2(x)%dat
    enddo

  else

    do x=xbeg,xend
     open(file_id_2_exp,file=trim(adjustl(pathname))//('Expt_tables/' &
                //trim(adjustl(S95_t2(x)%expt))//'tables/' &
                //trim(adjustl(S95_t2(x)%expt))//'tables2/' &
                //trim(adjustl(filename(x)))//'_pred.txt'))            
     open(file_id_2_obs,file=trim(adjustl(pathname))//('Expt_tables/' &
                //trim(adjustl(S95_t2(x)%expt))//'tables/' &
                //trim(adjustl(S95_t2(x)%expt))//'tables2/' &
                //trim(adjustl(filename(x)))//'_obs.txt'))                                 

     ! fill S95 from file
     ! row 0 and column 0 in LEP file contain higgs masses 
     ! and (0,0) ie top left set to -100
     ! so avoid them         
     allocate(testrow(0:S95_t2(x)%nx1))

     do k=lbound(file_id_arr,dim=1),ubound(file_id_arr,dim=1)
       read(file_id_arr(k),*)( testrow(i), i=0,S95_t2(x)%nx1 )
       if((testrow(0)+100.0D0).gt.small)stop'error in initializetables2 (d)' !top left number should be -100
       do i=1,S95_t2(x)%nx1      
         if( abs(testrow(i)- (S95_t2(x)%xmin1 + dble(i-1)*S95_t2(x)%sep1) ).gt.small*S95_t2(x)%sep1 )then
            write(*,*)S95_t2(x)%id,testrow(i),(S95_t2(x)%xmin1 + dble(i-1)*S95_t2(x)%sep1)
            stop'error in initializetables2 (e)'  
         endif
       enddo
     enddo
     deallocate(testrow)
  
     do j=1,S95_t2(x)%nx2  
      read(file_id_2_exp,*) dummy, ( S95_t2(x)%dat(j,i,2), i=1,S95_t2(x)%nx1 )
      if( abs(dummy- (S95_t2(x)%xmin2 + dble(j-1)*S95_t2(x)%sep2) ).gt.small*S95_t2(x)%sep2 ) then
       stop'error in initializetables2 (f)'
      endif
      read(file_id_2_obs,*) dummy, ( S95_t2(x)%dat(j,i,1), i=1,S95_t2(x)%nx1 )   
      if( abs(dummy- (S95_t2(x)%xmin2 + dble(j-1)*S95_t2(x)%sep2) ).gt.small*S95_t2(x)%sep2 ) then
       stop'error in initializetables2 (g)'                 
      endif
     end do
   
     maxdatval=S95_t2(x)%maxdatval
     if( maxdatval .gt. 0.0D0 )then
       ! set entries .ge. S95_t2(x)%maxdatval to (-4): they will not be relevent
       where(  S95_t2(x)%dat .ge. maxdatval )  S95_t2(x)%dat= - 4.0D0                        
     endif

     close(file_id_2_exp)
     close(file_id_2_obs)       
    enddo   

    rewind(file_id_common2)
#ifndef WEBVERSION
    do x=xbeg,xend
     write(file_id_common2)S95_t2(x)%dat
    enddo
#endif

  endif

  close(file_id_common2)
  deallocate(filename)      
                    
 end subroutine initializetables2
 !***********************************************************
 function t2elementnumberfromid(t2,id)
  !--------------------------------------input
  type(table2), intent(in) :: t2(:)
  integer, intent(in) :: id
  !-----------------------------------function
  integer :: t2elementnumberfromid
  !-----------------------------------internal
  integer :: n,x
  !-------------------------------------------

  n=0
  do x=lbound(t2,dim=1),ubound(t2,dim=1)
   if(t2(x)%id.eq.id)then
    n=n+1
    t2elementnumberfromid=x
   endif
  enddo

  if(n.ne.1)stop'problem in function t2elementnumberfromid 1'

 end function t2elementnumberfromid
 !*********************************************************** 
 subroutine fill_slices_t1_from_slices_of_t2(t2,v1orv2,xy_selection,ftype_selection,slices_t1)
 ! if this subroutine is used,  
 ! don't forget to deallocate slices_t1(x)%dat at some point
 !*********************************************************** 
 implicit none 
  !--------------------------------------input
  type(table2), intent(in) :: t2
  integer, intent(in) :: v1orv2
  integer, intent(in) :: xy_selection(:)
  integer, intent(in) :: ftype_selection(:)
  !-------------------------------------output
  type(table1) :: slices_t1(:)  !i.e. 2 slices
  !-----------------------------------internal
  integer :: i,j,k,n
  integer :: n_ftype_selection
  !-------------------------------------------
  n_ftype_selection=ubound(ftype_selection,dim=1)

  do n=lbound(ftype_selection,dim=1),n_ftype_selection
     if(ftype_selection(n).lt.lbound(t2%dat,dim=3))stop'problem in fill_slices_t1_from_slices_of_t2 3a'
     if(ftype_selection(n).gt.ubound(t2%dat,dim=3))stop'problem in fill_slices_t1_from_slices_of_t2 3b'
  enddo

  if(lbound(xy_selection,dim=1).ne.lbound(slices_t1,dim=1))then
     stop'problem in fill_slices_t1_from_slices_of_t2 1a'
  endif
  if(ubound(xy_selection,dim=1).ne.ubound(slices_t1,dim=1))then
     stop'problem in fill_slices_t1_from_slices_of_t2 1b'
  endif

  select case(v1orv2)
  case(1)

    do n=lbound(slices_t1,dim=1),ubound(slices_t1,dim=1)

     if(xy_selection(n).lt.lbound(t2%dat,dim=1))stop'problem in fill_slices_t1_from_slices_of_t2 4a'
     if(xy_selection(n).gt.ubound(t2%dat,dim=1))stop'problem in fill_slices_t1_from_slices_of_t2 4b'

     slices_t1(n)%id          =  t2%id    
     slices_t1(n)%nx          =  t2%nx1  
     slices_t1(n)%xmax        =  t2%xmax1 
     slices_t1(n)%xmin        =  t2%xmin1
     slices_t1(n)%sep         =  t2%sep1 
     slices_t1(n)%deltax      =  t2%deltax
  
     allocate( slices_t1(n)%dat(slices_t1(n)%nx,n_ftype_selection) ) 
     slices_t1(n)%dat = -1.0D0

     do i=1,slices_t1(n)%nx
       do k=1,n_ftype_selection
        slices_t1(n)%dat(i,k)=t2%dat(xy_selection(n),i,ftype_selection(k))
       enddo
     enddo

    enddo
  case(2)

    do n=lbound(slices_t1,dim=1),ubound(slices_t1,dim=1)

     if(xy_selection(n).lt.lbound(t2%dat,dim=2))stop'problem in fill_slices_t1_from_slices_of_t2 4aa'
     if(xy_selection(n).gt.ubound(t2%dat,dim=2))stop'problem in fill_slices_t1_from_slices_of_t2 4bb'

     slices_t1(n)%id          =  t2%id    
     slices_t1(n)%nx          =  t2%nx2  
     slices_t1(n)%xmax        =  t2%xmax2 
     slices_t1(n)%xmin        =  t2%xmin2
     slices_t1(n)%sep         =  t2%sep2 
     slices_t1(n)%deltax      =  t2%deltax
  
     allocate( slices_t1(n)%dat(slices_t1(n)%nx,n_ftype_selection) ) 
     slices_t1(n)%dat = -1.0D0

     do j=1,slices_t1(n)%nx
       do k=1,n_ftype_selection
        slices_t1(n)%dat(j,k)=t2%dat(j,xy_selection(n),ftype_selection(k))
       enddo
     enddo

    enddo
  case default
   stop'problem in fill_slices_t1_from_slices_of_t2 5'
  end select

 end subroutine fill_slices_t1_from_slices_of_t2
 !***********************************************************

 !*********************************************************** 
 subroutine fill_t1_from_t2(t2,v1orv2,xy_selection,ftype_selection,t1)
 ! if this subroutine is used,  
 ! don't forget to deallocate slices_t1(x)%dat at some point
 !*********************************************************** 
 implicit none 
  !--------------------------------------input
  type(table2), intent(in) :: t2
  integer, intent(in) :: v1orv2
  integer, intent(in) :: xy_selection
  integer, intent(in) :: ftype_selection(:)
  !-------------------------------------output
  type(table1) :: t1  
  !-----------------------------------internal
  integer :: i,j,k,n
  integer :: n_ftype_selection
  !-------------------------------------------
  n_ftype_selection=ubound(ftype_selection,dim=1)

  do n=lbound(ftype_selection,dim=1),n_ftype_selection
     if(ftype_selection(n).lt.lbound(t2%dat,dim=3))stop'problem in fill_t1_from_t2 3a'
     if(ftype_selection(n).gt.ubound(t2%dat,dim=3))stop'problem in fill_t1_from_t2 3b'
  enddo

  t1%id          =  t2%id 
  t1%deltax      =  t2%deltax

  select case(v1orv2)
  case(1)

     if(xy_selection.lt.lbound(t2%dat,dim=1))stop'problem in fill_t1_from_t2 4a'
     if(xy_selection.gt.ubound(t2%dat,dim=1))stop'problem in fill_t1_from_t2 4b'
  
     t1%nx          =  t2%nx1  
     t1%xmax        =  t2%xmax1 
     t1%xmin        =  t2%xmin1
     t1%sep         =  t2%sep1 
  
     allocate( t1%dat(t1%nx,n_ftype_selection) ) 
     t1%dat = -1.0D0

     do i=1,t1%nx
       do k=1,n_ftype_selection
        t1%dat(i,k)=t2%dat(xy_selection,i,ftype_selection(k))
       enddo
     enddo

  case(2)

     if(xy_selection.lt.lbound(t2%dat,dim=2))stop'problem in fill_t1_from_t2 4aa'
     if(xy_selection.gt.ubound(t2%dat,dim=2))stop'problem in fill_t1_from_t2 4bb'

     t1%nx          =  t2%nx2  
     t1%xmax        =  t2%xmax2 
     t1%xmin        =  t2%xmin2
     t1%sep         =  t2%sep2 
  
     allocate( t1%dat(t1%nx,n_ftype_selection) ) 
     t1%dat = -1.0D0

     do j=1,t1%nx
       do k=1,n_ftype_selection
        t1%dat(j,k)=t2%dat(j,xy_selection,ftype_selection(k))
       enddo
     enddo

  case default
   stop'problem in fill_t1_from_t2 5'
  end select

 end subroutine fill_t1_from_t2
 !*********************************************************** 
end module S95tables_type2
!************************************************************
