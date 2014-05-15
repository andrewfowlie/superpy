! This file is part of HiggsBounds
!  -KW
!****************************************************************** 
module interpolate
!******************************************************************
 implicit none
  integer,parameter :: n_points_max=50 ! for subroutine interpolate_1D_inv
 
 contains

 !******************************************************************       
 subroutine interpolate_tabletype1(x,t1,datcomp,interpol)
 !******************************************************************
 ! finds the value of 'interpol' at x by interpolating between points in table t1
 ! t1 is associated with a process involving only one Higgs
 ! datcomp indicates whether the predicted or observed experimental result is used
 ! all valid data in table must be >= 0 as
 ! negative numbers label invalid points
 !****************************************************************** 
  use S95tables_type1
  implicit none
  !--------------------------------------input
  double precision, intent(in) :: x
  integer, intent(in) :: datcomp
  type(table1), intent(in) :: t1
  !-------------------------------------output
  double precision :: interpol
  !-----------------------------------internal
  double precision :: c_x(2),c_f(2)
  logical :: intable
  !-------------------------------------------        

  call findcorners_tabletype1(x,t1,datcomp,intable,c_x,c_f)      

  if(intable)then
   interpol=interpol1D(x, &
                    &  c_x(1),c_f(1),  &
                    &  c_x(2),c_f(2) )
  else                 
   interpol=-1.0D0                
  endif
                
 end subroutine interpolate_tabletype1
 
 !******************************************************************       
 subroutine findcorners_tabletype1(x,t1,datcomp,intable,c_x,c_f)
 !******************************************************************
 ! finds relevant values in t1 for finding f(x) with interpolation 
 !
 ! first takes a table of type 1 (t1)
 ! then works out whether x is in the x-range of the table. 
 ! If not, sets intable=.False.
 ! If so, then finds the x-values above and below x, calls them c_x(1),c_x(2)
 ! and finds the corresponding f-values c_f(1),c_f(2)
 ! negative values of f(x) label invalid data points. If an invalid point 
 ! is needed for the interpolation, then intable is set to .False.
 !****************************************************************** 
  use usefulbits, only :small
  use S95tables_type1
  implicit none
  !--------------------------------------input
  double precision, intent(in) :: x
  integer, intent(in) :: datcomp
  type(table1), intent(in) :: t1
  !-------------------------------------output
  logical :: intable
  double precision :: c_x(2),c_f(2)
  !-----------------------------------internal
  integer :: ilow
  double precision :: x_below
  double precision :: xmin,xmax,sep
  !-------------------------------------------        

  intable=.True.                
  c_x=-1.0D0
  c_f=-1.0D0
                                                                                                   
  xmin=t1%xmin
  xmax=t1%xmax
                 
  ! check if mass is within mass range of table:
  if(.not. ( (x .ge. xmin-small).and.(x .le. xmax+small) )  )then ! written in convoluted way to get the NaNs
   intable=.False.
  else                
   sep=t1%sep                 
                    
   ilow=int((x-xmin)/sep)+1
                 
   x_below=dble(ilow-1)*sep+xmin

   if(abs(x_below-x).lt.small)then !x is the same as x_below 
    c_x(1)=x_below
    c_x(2)=c_x(1)
    c_f(1)=t1%dat(ilow,datcomp)
    c_f(2)=c_f(1)
   else
    c_x(1)=x_below
    c_x(2)=x_below+sep
    c_f(1)=t1%dat(ilow,datcomp)
    c_f(2)=t1%dat(ilow+1,datcomp)
   endif
                     
   if((c_f(1).lt.0.0D0).or.(c_f(2).lt.0.0D0))then
    intable=.False.                 
   endif
                                                                   
  endif 

                
 end subroutine findcorners_tabletype1

 !******************************************************************
 subroutine interpolate_tabletype2(x,y,t2,datcomp,interpol)
 !******************************************************************
 ! finds the value of 'interpol' by interpolating between points in table t2
 ! t2 is associated with a process involving two Higgs
 ! datcomp indicates whether the predicted or observed experimental result is used
 ! negative values of f(x) label invalid data points. If an invalid point 
 ! is needed for the interpolation, then intable is set to .False.
 !****************************************************************** 
  use S95tables_type2          
  implicit none
  !--------------------------------------input
  double precision, intent(in) :: y,x
  integer, intent(in) :: datcomp
  type(table2), intent(in) :: t2
  !-------------------------------------output
  double precision :: interpol
  !-----------------------------------internal
  logical :: intable
  double precision :: c_x(3),c_y(3),c_f(3)
  !-------------------------------------------              

  call findcorners_tabletype2(x,y,t2,datcomp,intable,c_x,c_y,c_f)

  if(intable)then
     interpol=interpol2D(x,y, &
                    &    c_x(1),c_y(1),c_f(1),&
                    &    c_x(2),c_y(2),c_f(2),& 
                    &    c_x(3),c_y(3),c_f(3))
  else
     interpol=-1.0D0                
  endif 

 end subroutine interpolate_tabletype2
 !*******************************************************                

 !******************************************************************
 subroutine findcorners_tabletype2(x,y,t2,datcomp,intable,c_x,c_y,c_f)
 !******************************************************************
 ! finds relevant values in t2 for finding f(x,y) with interpolation 
 !
 ! first takes a table of type 2 (t2)
 ! then works out whether x,y is in the x-range and y-range of the table. 
 ! ...if not, sets intable=.False.
 ! ...if so, then finds the x-values, and y-values around (x,y), 
 !          calls them c_x(1),c_x(2),c_x(3),c_y(1),c_y(2),c_y(3)
 ! and finds the corresponding f(x,y): c_f(1),c_f(2),c_f(3)
 ! Then checks that all f are positive (negative values in tabletype2 lable invalid data points)
 ! ...if not, sets intable=.False.
 !****************************************************************** 
  use S95tables_type2                
  use usefulbits, only :small
  implicit none
  !--------------------------------------input
  double precision, intent(in) :: y,x
  integer, intent(in) :: datcomp
  type(table2), intent(in) :: t2
  !-------------------------------------output
  logical :: intable
  double precision :: c_x(3),c_y(3),c_f(3)
  !-----------------------------------internal
  integer :: ilow,jlow,n
  integer :: ny,nx
  double precision :: ymin,ymax,xmin,xmax,sepy,sepx
  double precision :: yhigh,ylow,xhigh,xlow
  double precision :: testtri
  double precision :: ybit,xbit
  integer :: c_i(3),c_j(3) ! c are the relevant corners of the interpolation region
  !-------------------------------------------
  
  intable=.True.                
  c_x=-1.0D0
  c_y=-1.0D0
  c_f=-1.0D0   

  ymin=t2%xmin2
  ymax=t2%xmax2
  xmin=t2%xmin1
  xmax=t2%xmax1                                

  ! check if mass is within mass range of table:

  if(    .not.((x.gt.-1.0D0).or.(x.le.0.0D0)))then! written in convoluted way to get the NaNs
    intable=.False.
  elseif(.not.((y.gt.-1.0D0).or.(y.le.0.0D0)))then! written in convoluted way to get the NaNs
    intable=.False.

  elseif( (x.lt.xmin).and.(abs(x-xmin)**2.0D0.gt.small) )then !x is too low to be in table
    intable=.False.
  elseif( (y.lt.ymin).and.(abs(y-ymin)**2.0D0.gt.small) )then !y is too low to be in table
    intable=.False.

  elseif( (x.gt.xmax).and.(abs(x-xmax)**2.0D0.gt.small) )then !x is too high to be in table
    intable=.False.
  elseif( (y.gt.ymax).and.(abs(y-ymax)**2.0D0.gt.small) )then !y is too high to be in table
    intable=.False.
                               
  else

   ny=t2%nx2
   nx=t2%nx1
   sepy=t2%sep2
   sepx=t2%sep1                                
 
  ! Interpolation: 
  ! points in table are treated as corners of flat triangular surfaces                 

   jlow=1+int( (y-ymin) /sepy )
   !ybit=((y-ymin)-dble(int(y-ymin)))/sepy
   ybit=(y-(ymin+dble(jlow-1)*sepy))/sepy
                 
   ilow=1+int( (x-xmin) /sepx )
   !xbit=((x-xmin)-dble(int(x-xmin)))/sepx                
   xbit=(x-(xmin+dble(ilow-1)*sepx))/sepx

   if((ybit*sepy)**2.0D0+(xbit*sepx)**2.0D0 .le. small)then! exactly on datapoint
     c_i(1)=ilow  ; c_j(1)=jlow  !i.e only need to use ilow,jlow
     c_i(2)=ilow  ; c_j(2)=jlow     
     c_i(3)=ilow  ; c_j(3)=jlow     
    
   else

     testtri=ybit/sqrt(ybit**2.0D0+xbit**2.0D0) 
     xlow =dble(ilow-1  )*sepx+xmin
     xhigh=dble(ilow-1+1)*sepx+xmin
     ylow =dble(jlow-1  )*sepy+ymin
     yhigh=dble(jlow-1+1)*sepy+ymin

     if((xbit*sepx)**2.0D0 .le. small)then !want linear interpolation along y axis
      c_i(1)=ilow  ; c_j(1)=jlow !i.e only need to use (ilow,jlow) and (ilow,jlow+1)
      c_i(2)=ilow  ; c_j(2)=jlow+1
      c_i(3)=ilow  ; c_j(3)=jlow+1

    !now sort out edge effects for tables where t2%needs_M2_gt_2M1=True
    !            x,i  -->                  
    !            o       o 
    !            *-    
    !            *--               upper square
    !            *---                  
    !            o----   o  
    !        |   * *--- B   
    !        v   * A *--           lower square
    !       y,j  *     *-    
    !            o*******o

     elseif(  t2%needs_M2_gt_2M1  &
       &   .and. &
       & (y.ge.2.0D0*x) &!y is above or on the line y=2x
       &   .and.&
       & ((yhigh+small).lt.2.0D0*xhigh) &!upper square
       &   .and. &
       & ((jlow+2).le.ny) &!all points of triangle for interpolation must be in table
       )then  
      c_i(1)=ilow    ; c_j(1)=jlow
      c_i(2)=ilow    ; c_j(2)=jlow+1
      c_i(3)=ilow+1  ; c_j(3)=jlow+2 

     elseif(  t2%needs_M2_gt_2M1 &
       &   .and. &
       & (y.ge.2.0D0*x) & !y is above or on the line y=2x
       &   .and. &
       & ((ylow+small).lt.2.0D0*xhigh).and.(testtri.le.0.7070D0) &!lower square, triangle B
       &   .and. &
       & ((jlow-1).gt.0) & !all points of triangle for interpolation must be in table
       )then 
      c_i(1)=ilow    ; c_j(1)=jlow-1 
      c_i(2)=ilow    ; c_j(2)=jlow
      c_i(3)=ilow+1  ; c_j(3)=jlow+1  

     elseif((ybit*sepy)**2.0D0 .le. small)then !want linear interpolation along x axis
      c_i(1)=ilow  ; c_j(1)=jlow  !i.e only need to use (ilow,jlow) and (ilow+1,jlow)
      c_i(2)=ilow+1; c_j(2)=jlow
      c_i(3)=ilow+1; c_j(3)=jlow

     else   
    !            x,i  -->                  
    !            o*******o  testtri will find out which triangle required point r is in
    !        |   * *   B *  i.e. triangle of orientation 'A' or 'B'  
    !        v   * A *   *
    !       y,j  *     * *    
    !            o*******o  
    !
    !            xbit
    !            <->                      
    !        ^   o*******o jlow,ilow+1  
    !    ybit|   * *     *    
    !        v   * r *   *
    !            *     * *    
    ! jlow+1,ilowo*******oilow+1,ilow+1   

      if(testtri.gt.0.7070D0)then !triangle orientation A    
     !1/sqrt(2)=0.707107:slight bias towards orientation A, good for diagonal edge of graph

      ! find 'height' (i.e. S95 value) of each corner of triangle
                   
      !         c3*       
      !           * *          
      !           * r * 
      !           *     *      
      !         c1*********c2

        c_j(1)=jlow+1 ;c_i(1)=ilow
        c_j(2)=jlow+1 ;c_i(2)=ilow+1
        c_j(3)=jlow   ;c_i(3)=ilow
                                                        
      else !triangle orientation B                
      !       c2*********c1 
      !           *     *   
      !             *   *  
      !               * *    
      !                 *c3  

        c_j(1)=jlow  ; c_i(1)=ilow+1
        c_j(2)=jlow  ; c_i(2)=ilow
        c_j(3)=jlow+1; c_i(3)=ilow+1

      endif
     endif
   endif

   c_x(:)=dble(c_i(:)  -1)*sepx+xmin  !x value of each corner
   c_y(:)=dble(c_j(:)  -1)*sepy+ymin                 !y value of each corner

   do n=1,ubound(c_f,dim=1)
     c_f(n)=t2%dat(c_j(n),c_i(n),datcomp) !S95 value of each corner
   enddo

   if( (c_f(1).lt.0.0D0).or.(c_f(2).lt.0.0D0).or.(c_f(3).lt.0.0D0) )then !check that all of the point needed for interpolation
    intable=.False.
   endif

  endif 

 end subroutine findcorners_tabletype2
 !*******************************************************     
 subroutine interpolate_tabletype3(x,y,z,t3,datcomp,interpol)
 ! interpolates in table t3 to find interpol
 ! are given by interpol=f(x,y,z)
 ! where f is stored in t3%dat(:,:,datcomp)
 !******************************************************************
  use S95tables_type2                
  use S95tables_type3            
  implicit none
  !--------------------------------------input
  double precision, intent(in) :: y,x,z
  integer, intent(in) :: datcomp
  type(table3), intent(in) :: t3
  !-------------------------------------output
  double precision :: interpol
  !-----------------------------------internal
  type(table2) :: slices_t2(2)
  integer :: a
  integer :: ftype_selection(1) !ftype_selection(1)=datcomp
  logical :: filledslice
  double precision :: interpol_array(1)
  double precision :: c_x(3,2) !x-values of 3 points where the value of f is known. 2nd element: above,below
  double precision :: c_y(3,2) !y-values of 3 points where the value of f is known. 2nd element: above,below
  double precision :: c_z(2) !z-values of 2 points where the value of f is known, both with same x-value and y-values
  double precision :: c_f(3,2)  !f(n,i) is f at (c_x(n,1),c_y(n,1),c_z(i))=(c_x(n,2),c_y(n,2),c_z(i)). 2nd element: above,below
  !-------------------------------------------
  ftype_selection(1)=datcomp

  call interpolate_tabletype3_longer(x,y,z,t3,ftype_selection,interpol_array, & 
         & slices_t2,filledslice,c_x,c_y,c_z,c_f)
  interpol=interpol_array(1)

  if(filledslice)then
   do a=1,2  
     deallocate(slices_t2(a)%dat)
   enddo
  endif

 end subroutine interpolate_tabletype3               
 !******************************************************************       
 subroutine interpolate_tabletype3_longer(x,y,z,t3,ftype_selection,interpol_array, & 
         & slices_t2,filledslice,c_x,c_y,c_z,c_f)
 ! interpolates in table t3 to find interpol
 ! are given by interpol=f(x,y,z)
 ! where f is stored in t3%dat(:,:,ftype_selection(:))
 ! if filledslice=.True., must remember to deallocate slices_t2(a)%dat
 ! c_f corresponds to ftype_selection(1) only
 !******************************************************************
  use S95tables_type2                
  use S95tables_type3                
  use usefulbits, only :small
  implicit none
  !--------------------------------------input
  double precision, intent(in) :: y,x,z
  type(table3), intent(in) :: t3
  integer, intent(in) :: ftype_selection(:) 
  !-------------------------------------output
  double precision :: interpol_array(:)
  type(table2) :: slices_t2(2)
  logical :: filledslice
  double precision :: c_x(3,2) !x-values of 3 points where the value of f is known. 2nd element: above,below
  double precision :: c_y(3,2) !y-values of 3 points where the value of f is known. 2nd element: above,below
  double precision :: c_z(2) !z-values of 2 points where the value of f is known, both with same x-value and y-values
  double precision :: c_f(3,2)  !f(n,i) is f at (c_x(n,1),c_y(n,1),c_z(i))=(c_x(n,2),c_y(n,2),c_z(i)). 2nd element: above,below
  !-----------------------------------internal
  double precision :: z_below,z_above
  integer :: ilow,a
  integer :: c_zi(2)
  !-------------------------------------------

  if(lbound(interpol_array, dim=1).ne.lbound(ftype_selection, dim=1))then
    stop'problem in interpolate_tabletype3_longer (a)'
  elseif(ubound(interpol_array, dim=1).ne.lbound(ftype_selection, dim=1))then
    stop'problem in interpolate_tabletype3_longer (b)'  
  endif

  ! check if mass is within z range of table:
  if(    .not. ( (z .ge. t3%zmin-small).and.(z .le. t3%zmax+small) )  )then !#1! written in convoluted way to get the NaNs
   interpol_array=-1.0D0      
   filledslice=.False.
   c_x=-1.0D0
   c_y=-1.0D0
   c_z=-1.0D0
   c_f=-1.0D0
   
  else              !#1
                    
   ilow=int((z-t3%zmin)/t3%zsep)+1
   z_below=dble(ilow-1)*t3%zsep+t3%zmin
   z_above=z_below+t3%zsep

   if(abs(z_below-z).lt.small)then !z is the same as z_below 
    c_zi= ilow
   elseif(abs(z_above-z).lt.small)then !z is the same as z_above                 
    c_zi= ilow+1
   else
    c_zi(1)= ilow
    c_zi(2)= ilow+1
   endif

   call fill_slices_t2_from_slices_of_t3(t3,c_zi,ftype_selection,slices_t2)
   
   a=1
   call interpolate_slices_t2_longer(x,y,z,slices_t2,a,interpol_array(a),c_x,c_y,c_z,c_f)

   do a=1+lbound(ftype_selection, dim=1),ubound(ftype_selection, dim=1)
     call interpolate_slices_t2(x,y,z,slices_t2,a,interpol_array(a))
   enddo

   filledslice=.True.  
     
  endif !#1

 end subroutine interpolate_tabletype3_longer
 !******************************************************************       
 subroutine interpolate_slices_t2(x,y,z,slices_t2,datcomp,interpol)
 ! interpolates in slices_t2 to find interpol
 ! are given by interpol=f(x,y,z)
 ! where f is stored in slices_t2%dat(:,:,datcomp)
 !******************************************************************
  use S95tables_type2                
  implicit none
  !--------------------------------------input
  double precision, intent(in) :: y,x,z
  integer, intent(in) :: datcomp
  type(table2), intent(in) :: slices_t2(:)
  !-------------------------------------output
  double precision :: interpol
  !-----------------------------------internal
  double precision :: c_x(3,2) !x-values of 3 points where the value of f is known. 2nd element: above,below
  double precision :: c_y(3,2) !y-values of 3 points where the value of f is known. 2nd element: above,below
  double precision :: c_z(2) !z-values of 2 points where the value of f is known, both with same x-value and y-values
  double precision :: c_f(3,2)  !f(n,i) is f at (c_x(n,1),c_y(n,1),c_z(i))=(c_x(n,2),c_y(n,2),c_z(i)). 2nd element: above,below
  !-------------------------------------------

  call interpolate_slices_t2_longer(x,y,z,slices_t2,datcomp,interpol,c_x,c_y,c_z,c_f)

 end subroutine interpolate_slices_t2
 !******************************************************************       
 subroutine interpolate_slices_t2_longer(x,y,z,slices_t2,datcomp,interpol,c_x,c_y,c_z,c_f)
 ! interpolates in slices_t2 to find interpol
 ! are given by interpol=f(x,y,z)
 ! where f is stored in slices_t2%dat(:,:,datcomp)
 !******************************************************************
  use S95tables_type2                
  use S95tables_type3                
  use usefulbits, only :small
  implicit none
  !--------------------------------------input
  double precision, intent(in) :: y,x,z
  integer, intent(in) :: datcomp
  type(table2), intent(in) :: slices_t2(:)
  !-------------------------------------output
  double precision :: interpol
  double precision :: c_x(3,2) !x-values of 3 points where the value of f is known. 2nd element: above,below
  double precision :: c_y(3,2) !y-values of 3 points where the value of f is known. 2nd element: above,below
  double precision :: c_z(2) !z-values of 2 points where the value of f is known, both with same x-value and y-values
  double precision :: c_f(3,2)  !f(n,i) is f at (c_x(n,1),c_y(n,1),c_z(i))=(c_x(n,2),c_y(n,2),c_z(i)). 2nd element: above,below
  !-----------------------------------internal
  logical :: intable,intableslice(2)
  double precision :: z_below,z_above
  integer :: i,ilow,a
  integer :: one_good_edge
  logical :: use_edge(3)
  logical :: pointisinsidetable
  !-------------------------------------------
  
  intable=.True.                
  intableslice=.True.
  c_x   =  -1.0D6                
  c_y   =  -1.0D6
  c_f   =  -1.0D6

  if(size(slices_t2,dim=1).ne.2)stop'wrong input to slices_t2 (1)'

  ! check if mass is within z range of table:

  if(    .not. ( (z .ge. slices_t2(1)%z-small).and.(z .le. slices_t2(2)%z+small) )  )then !#1! written in convoluted way to get the NaNs
   intable=.False.                 
  else                !#1
                    
   ilow=int((z-slices_t2(1)%z)/(slices_t2(2)%z-slices_t2(1)%z))+1
   z_below=dble(ilow-1)*(slices_t2(2)%z-slices_t2(1)%z)+slices_t2(1)%z
   z_above=z_below+(slices_t2(2)%z-slices_t2(1)%z)

   if(abs(z_below-z).lt.small)then !z is the same as z_below 
    c_z=z_below
    intableslice(2)=.False.
   elseif(abs(z_above-z).lt.small)then !z is the same as z_above 
    c_z=z_above
    intableslice(1)=.False.
   else
    c_z(1)=z_below
    c_z(2)=z_above
   endif

   pointisinsidetable=.True.
   do a=1,2  
    if(intableslice(a))then
       call findcorners_tabletype2(x,y,slices_t2(a),datcomp,pointisinsidetable,c_x(:,a),c_y(:,a),c_f(:,a))
    endif
   enddo

   if((.not.intableslice(1)).and.(.not.intableslice(2)))then !#3
    intable=.False.
   elseif(.not.pointisinsidetable)then!#3
   else!#3
     if( .not.intableslice(1))then !if one of the slices are not needed, set values at corners to be the same as the slice that is needed
      call copyslice(1,2)
     elseif( .not.intableslice(2))then
      call copyslice(2,1)
     endif

     one_good_edge=-1 ! edges are constants in (x,y)
     use_edge=.True.
     do i=1,3
       if((abs(c_x(i,1)-c_x(i,2)).gt.small).and.(abs(c_y(i,1)-c_y(i,2)).gt.small) )then ! (c_x(1,i),c_y(1,i)) is not the same point as (c_x(2,i),c_y(2,i))
         use_edge(i) = .False.
       else
         one_good_edge=i
       endif
     enddo

     if(one_good_edge.lt.0)then !#4! no good edges. need at least one for interpolation
                                ! should only occur when there's only one valid point on each slice and they don't line up
      intable=.False.
     else!#4
     

        if(use_edge(1).and.use_edge(2).and.use_edge(3))then!#5 ! all three edges can be used

        else!#5
         do i=1,3 ! if one of the edges can't be used, need to replace it with a good edge
           if(.not.use_edge(i))then
            call copyedge(i,one_good_edge)
           endif
         enddo
        endif!#5
      endif!#4
     endif!#3
  endif  !#1

  !write(*,*)'hello cx,cy,cz,cf',c_x,c_y,c_z,c_f
  if(intable.and.(minval(c_f).ge.0.0D0))then !#2

    interpol= interpol3D( x,y,z, &
                    &    c_x(1,1),c_y(1,1),c_z,c_f(1,:),&
                    &    c_x(2,1),c_y(2,1),c_z,c_f(2,:),& 
                    &    c_x(3,1),c_y(3,1),c_z,c_f(3,:))

  else                  !#2
   interpol=-1.0D0                
  endif !#2

  contains

  subroutine copyedge(aa,bb)
    integer,intent(in) :: aa,bb
      c_x(aa,:)   =  c_x(bb,:)                
      c_y(aa,:)   =  c_y(bb,:) 
      c_f(aa,:)   =  c_f(bb,:)
  end subroutine copyedge

  subroutine copyslice(aa,bb)
   integer,intent(in) :: aa,bb
      c_x(:,aa)   =  c_x(:,bb)                
      c_y(:,aa)   =  c_y(:,bb) 
      c_z(aa)     =  c_z(bb) 
      c_f(:,aa)   =  c_f(:,bb)
  end subroutine copyslice

 end subroutine interpolate_slices_t2_longer
 !*******************************************************
 function interpol1D(x,x1,f1,x2,f2)
 !*******************************************************
  use usefulbits, only :small
   !-----------------------------------function
   double precision interpol1D ! function f at x
   !--------------------------------------input
   double precision, intent(in) :: x          
   double precision, intent(in) :: x1,x2      !x-values where the value of f is known
   double precision, intent(in) :: f1         !f1 is f at x1
   double precision, intent(in) :: f2         !f2 is f at x2
   !-------------------------------------------

   if(abs(x-x1)**2.0 .lt. small)then
    interpol1D=f1
   elseif(abs(x-x2)**2.0 .lt. small)then
    interpol1D=f2
   elseif(abs(x2-x1)**2.0 .gt. small)then
    interpol1D=f1 +(f2-f1)/(x2-x1)*(x-x1) 
    if(  ((x1-x).lt.small) .eqv. ((x-x2).le.small) )then !shorter than  ( ((x1-x).lt.small) .and. ((x-x2).le.small) ) .or. ( ((x2-x).le.small) .and. ((x-x1).lt.small) ) 
     ! x is between x1 and x2
    else 
     ! x is not between x1 and x2 i.e. extrapolation, not interpolation
     write(*,*)'x,x1,x2',x,x1,x2
     write(*,*)'Warning: function interpol1D may not be reliable'
    endif
   else ! x1 and x2 are the same point, but x is somewhere else
     write(*,*)'x,x1,x2',x,x1,x2
     stop'error in function interpol1D' 
   endif

 end function interpol1D

 !*******************************************************                
 function interpol2D(x,y,x1,y1,f1,x2,y2,f2,x3,y3,f3)
  use usefulbits, only :small
  ! interpolation within a triangle on the xy plane 
  ! note that it is *not* necessary to have x1!=x2!=x3
  !                                     and y2!=y2!=y3
  ! f1=f(x1,y1)   f2=f(x2,y2)     
  !      *************
  !       *   P    *          P is (x,y), interpol2D is the value of f at P  
  !         *   *
  !           *f3=f(x3,y3)  

   !-----------------------------------function
   double precision interpol2D ! function f at x,y
   !--------------------------------------input
   double precision, intent(in) :: x,y          
   double precision, intent(in) :: x1,x2,x3,y1,y2,y3 !x-values and y-values of points where the value of f is known
   double precision, intent(in) :: f1         !f1 is f at (x1,y1)
   double precision, intent(in) :: f2         !f2 is f at (x2,y2)
   double precision, intent(in) :: f3         !f3 is f at (x3,y3)
   !-----------------------------------internal
   double precision denom   
   double precision r,r1,r2,r3   
   logical point_is_inside_Tri
   !-------------------------------------------                
   denom= (x2-x1)*(y3-y1)-(y2-y1)*(x3-x1)

   if(    ((x1-x )**2.0D0 + (y1-y )**2.0D0).lt. small )then ! (x,y) is at the same place as (x1,y1)
     interpol2D=f1    
     !write(*,*)'hello: case 1a'
   elseif(((x2-x )**2.0D0 + (y2-y )**2.0D0).lt. small )then ! (x,y) is at the same place as (x2,y2)
     interpol2D=f2    
     !write(*,*)'hello: case 1b'
   elseif(((x3-x )**2.0D0 + (y3-y )**2.0D0).lt. small )then ! (x,y) is at the same place as (x3,y3)
     interpol2D=f3  
     !write(*,*)'hello: case 1c'

   elseif(abs(denom).lt.small)then  
     if(  ((x1-x2)**2.0D0 + (y1-y2)**2.0D0 &
       & + (x1-x3)**2.0D0 + (y1-y3)**2.0D0).lt. small )then ! (x1,y1),(x2,y2),(x3,y3) are in the same place but (x,y) is at a different place
       stop'error in function interpol2D (case 1d)' 

     elseif( ((x1-x2)**2.0D0 + (y1-y2)**2.0D0) .lt. small )then ! (x1,y1) and (x2,y2) are the same points
       if( abs((x2-x3)*(y-y3)-(y2-y3)*(x-x3)) .lt.sqrt(small) .or. &
          & (  ((x2-x3)**2.0D0.lt.small) .and. ((x-x3)**2.0D0.lt.small) ).or. &
          & (  ((y2-y3)**2.0D0.lt.small) .and. ((y-y3)**2.0D0.lt.small) ) &
          &    )then !the points are all in a line
         !write(*,*)'hello: case 2ai'

         ! basis change 
         r=  sqrt( (x -x2)**2.0D0 + (y -y2)**2.0D0)
         r2= 0.0D0
         r3= sqrt( (x3-x2)**2.0D0 + (y3-y2)**2.0D0)

         interpol2D= interpol1D(r,r2,f2,r3,f3)
       else
         stop'error in function interpol2D (case 2aii)' 
       endif
     elseif( ((x1-x3)**2.0D0 + (y1-y3)**2.0D0) .lt. small )then ! (x1,y1) and (x3,y3) are the same points
       if( abs((x2-x1)*(y-y1)-(y2-y1)*(x-x1)) .lt.sqrt(small) .or. &
          & (  ((x2-x1)**2.0D0.lt.small) .and. ((x-x1)**2.0D0.lt.small) ).or. &
          & (  ((y2-y1)**2.0D0.lt.small) .and. ((y-y1)**2.0D0.lt.small) ) &
          &    )then !the points are all in a line
         !write(*,*)'hello: case 2bi'
         ! basis change 
         r=  sqrt( (x -x1)**2.0D0 + (y -y1)**2.0D0)
         r1= 0.0D0
         r2= sqrt( (x2-x1)**2.0D0 + (y2-y1)**2.0D0)

         interpol2D= interpol1D(r,r1,f1,r2,f2)
       else
         stop'error in function interpol2D (case 2bii)' 
       endif

     elseif( ((x2-x3)**2.0D0 + (y2-y3)**2.0D0) .lt. small )then ! (x2,y2) and (x3,y3) are the same points

       if(  (  abs((x2-x1)*(y-y1)-(y2-y1)*(x-x1)) .lt.sqrt(small) ) .or. &
          & (  ((x2-x1)**2.0D0.lt.small) .and. ((x-x1)**2.0D0.lt.small) ).or. &
          & (  ((y2-y1)**2.0D0.lt.small) .and. ((y-y1)**2.0D0.lt.small) ) &
          &    )then !the points are all in a line
         !write(*,*)'hello: case 2ci'
         ! basis change 
         r=  sqrt( (x -x1)**2.0D0 + (y -y1)**2.0D0)
         r1= 0.0D0
         r3= sqrt( (x3-x1)**2.0D0 + (y3-y1)**2.0D0)

         interpol2D= interpol1D(r,r1,f1,r3,f3)
       else
         stop'error in function interpol2D (case 2cii)'
       endif

     elseif( abs((x2-x1)*(y-y1)-(y2-y1)*(x-x1)) .lt.small )then !the points are all in a line, but none are in exactly the same place
         !first, need to work out which 2 points to use
         !then, use interpol1D
         stop'error: have not implemented this bit yet (case 2d)'
     else !(x1,y1),(x2,y2),(x3,y3) are in a line but (x,y) is not in that line
         stop'error in function interpol2D (case 2e)'
     endif
   else !case 3
   ! write(*,*)'hello: case 3'

   ! define a plane going through (x1,y1),(x2,y2),(x3,y3)
   ! then find 'height' above xy plane at point (x,y)

     interpol2D=  -(x-x1)*( (y2-y1)*(f3-f1) - (f2-f1)*(y3-y1) )/denom &
        &         +(y-y1)*( (x2-x1)*(f3-f1) - (f2-f1)*(x3-x1) )/denom &
        &         +f1

     ! check that (x,y) is in triangle with corners at (x1,y1),(x2,y2),(x3,y3)
     call check_if_point_is_inside_Tri
     if(point_is_inside_Tri.eqv..False.) write(*,*)'Warning: function interpol 2D may not be reliable' !i.e. extrapolation, not interpolation

   endif
   contains
   !----------------------------------
   function detof2vec(a,b)
    implicit none 
    double precision detof2vec
    double precision, intent(in) :: a(2),b(2)
   
    detof2vec=a(1)*b(2)-b(1)*a(2)

   end function detof2vec
   !----------------------------------
   subroutine check_if_point_is_inside_Tri
   !checks if point (x,y) is inside triangle (x1,y1),(x2,y2),(x3,y3)
   !n.b. have already checked that (x1,y1),(x2,y2),(x3,y3) forms a triangle rather than a point or a line
   !using formulae from http://mathworld.wolfram.com/TriangleInterior.html
    implicit none 
    double precision a,b,eps
    double precision v(2),v0(2),v1(2),v2(2)  
   
    v(1)=x
    v(2)=y

    v0(1)=x1
    v0(2)=y1

    v1(1)=x2-x1
    v1(2)=y2-y1

    v2(1)=x3-x1
    v2(2)=y3-y1

    a=   (detof2vec(v,v2)-detof2vec(v0,v2))/detof2vec(v1,v2)
    b= - (detof2vec(v,v1)-detof2vec(v0,v1))/detof2vec(v1,v2)

    ! (a+b).le.(1.0D0+eps) because want points near to diagonal edge as well
    eps=0.002 ! need to calculate eps properly from angle in subroutine findcorners_tabletype2
    if( (a+small.ge.0.0D0) .and. (b+small.ge.0.0D0) .and. ( (a+b).le.(1.0D0+eps) ) )then 
      point_is_inside_Tri=.True.
    else
      point_is_inside_Tri=.False.
    endif

   end subroutine check_if_point_is_inside_Tri
   !----------------------------------

  end function interpol2D

 !*******************************************************
 function interpol3D(x,y,z,x1,y1,z1,f1,x2,y2,z2,f2,x3,y3,z3,f3)
!  use usefulbits, only :small
  ! interpolation within a triangular-based prism, where 
  ! the z axis is parallel to the 'height'
  ! note that it is *not* necessary to have x1!=x2!=x3
  !                                     and y1!=y2!=y3
  !
  !  P is (x,y,z), interpol3D is the value of f at P  
  ! 
  ! f1(2)=f(x1,y1,z1(2))        f2(2)=f(x2,y2,z2(2))       z-axis
  !                * * * * * * * * * *                       ^
  !                |   *           * |                       |
  !                |       *     *   |                       |
  !                |           *f3(2)=f(x3,y3,z3(2))         |
  !                |           |     |
  !f1mid=f(x1,y1,z)|...........|.....|f2mid=f(x2,y2,z)
  !                |   . P     |   . |
  !                |       .   | .   |
  !                |           |f3mid=f(x3,y3,z)
  !                |           |     |
  ! f1(1)=         * * * * * * * * * * f2(1)=f(x2,y2,z2(1))
  !f(x1,y1,z1(1))      *       |   *
  !                        *   | *
  !                            *f3(1)=f(x3,y3,z3(1))
  !

   !-----------------------------------function
   double precision interpol3D ! function f at x,y,z
   !--------------------------------------input
   double precision, intent(in) :: x,y,z          
   double precision, intent(in) :: x1,x2,x3 !x-values of 2 points where the value of f is known
   double precision, intent(in) :: y1,y2,y3 !y-values of 2 points where the value of f is known
   double precision, intent(in) :: z1(2),z2(2),z3(2) !z-values of 2 points where the value of f is known, both with same x-value and y-value 
   ! f1 is declared as f1(:) rather than f1(2) because ifort complains otherwise (same with f2,f3)
   double precision, intent(in) :: f1(:)         !f1(i) is f at (x1,y1,z1(i))
   double precision, intent(in) :: f2(:)         !f2(i) is f at (x2,y2,z2(i))
   double precision, intent(in) :: f3(:)         !f3(i) is f at (x3,y3,z3(i))
   !-----------------------------------internal
   double precision f1mid,f2mid,f3mid
   !------------------------------------------- 

   if(    (lbound(f1,dim=1).ne.1) .or. (lbound(f2,dim=1).ne.1) .or. (lbound(f3,dim=1).ne.1) )then
     stop'problem in function interpol3D (1)'
   elseif((ubound(f1,dim=1).ne.2) .or. (ubound(f2,dim=1).ne.2) .or. (ubound(f3,dim=1).ne.2) )then  
     stop'problem in function interpol3D (2)'      
   endif

   f1mid = interpol1D(z,z1(1),f1(1),z1(2),f1(2))
   f2mid = interpol1D(z,z2(1),f2(1),z2(2),f2(2))   
   f3mid = interpol1D(z,z3(1),f3(1),z3(2),f3(2))

   interpol3D=interpol2D( x,y,          &
                    &     x1,y1,f1mid,  &
                    &     x2,y2,f2mid,  & 
                    &     x3,y3,f3mid   )

 end function interpol3D
 !******************************************************************       
 subroutine interpolate_1D_inv(targt,t1,datcomp,m_in,m_interpol,m_interpol_nearest,n_points)
 !******************************************************************
 ! finds the masses at which the experimental limit equals targt.
 ! m_interpol will be the closest of these masses to m_in.
 ! the m_in should be positive
 ! negative m_interpol labels point outside table
 ! treats negative values of t1%dat as invalid
 !need to tidy up a bit
 !****************************************************************** 
  use S95tables_type1
  use usefulbits, only :small
  implicit none
  !--------------------------------------input
  double precision, intent(in) :: targt
  integer, intent(in) :: datcomp
  type(table1), intent(in) :: t1
  double precision, intent(in) :: m_in
  !-------------------------------------output
  double precision, intent(out) :: m_interpol(:),m_interpol_nearest
  integer, intent(out) ::  n_points
  !-----------------------------------internal
  integer :: i,c
  double precision :: m_temp,m_i,m_iplus1,f_i,f_iplus1
  double precision,allocatable :: dat(:)
  !-------------------------------------------        

  if(size(m_interpol,dim=1).ne.n_points_max)stop'wrong input to interpolate_1D_inv 1'

 !table 1
 !fact=   1.2001890350471742     
 !M0=   115.18553525008736 
  m_interpol= -1.0D8
  m_interpol_nearest= -1.0D8

  allocate(dat(lbound(t1%dat,dim=1):ubound(t1%dat,dim=1)))
  dat=t1%dat(:,datcomp)

  c=0

  if(    minval(dat, mask=(dat.ge.0.0D0) ).gt.targt)then! off the upper end of the table
   m_interpol= - 3.0D0
  elseif(maxval(dat, mask=(dat.ge.0.0D0) ).lt.targt)then! off the lower end of the table
   m_interpol= - 4.0D0
  elseif(m_in.gt.t1%xmax)then !m_in should also be in table
   m_interpol= - 4.0D0
   stop'wrong input to interpolate_1D_inv a' 
  elseif(m_in.lt.t1%xmin)then
   m_interpol= - 4.0D0
   stop'wrong input to interpolate_1D_inv b'

  else
   do i=1,t1%nx-1 ! the case where targt is exactly on last data point is dealt with separately below

     f_i      = dat(i  )
     f_iplus1 = dat(i+1)     !notation: "iplus1"

     if(abs(f_i-targt).lt.small)then !targt is exactly on data point
       m_temp=dble(i-1)*t1%sep + t1%xmin

       c=c+1;if(c.le.n_points_max)m_interpol(c)=m_temp

     elseif( min(f_i,f_iplus1) .lt. 0.0D0 )then
       m_temp=-1.0D8

     elseif(  ( targt.lt. max( f_i,f_iplus1 ) ) .and. &
          &   ( targt.gt. min( f_i,f_iplus1 ) ) )then

        m_i      = dble(i  -1)*t1%sep + t1%xmin 
        m_iplus1 = dble(i+1-1)*t1%sep + t1%xmin

        m_temp=interpol1D(targt,f_i,m_i,f_iplus1,m_iplus1)

        c=c+1;if(c.le.n_points_max)m_interpol(c)=m_temp

     else
       m_temp=-1.0D8
     endif

     if((abs(m_temp-m_in).lt.abs(m_interpol_nearest-m_in)).and.(m_temp.ge.0.0D0))then
       m_interpol_nearest=m_temp   !want m_intepol to be as near to m_in as possible
     endif

   enddo

   if(abs( dat(t1%nx) - targt ).lt.small)then !targt is exactly on last data point
     m_temp=dble(t1%nx-1)*t1%sep + t1%xmin

     c=c+1;if(c.le.n_points_max)m_interpol(c)=m_temp

     if((abs(m_temp-m_in).lt.abs(m_interpol_nearest-m_in)).and.(m_temp.ge.0.0D0))then
       m_interpol_nearest=m_temp   !want m_intepol to be as near to m_in as possible
     endif
   endif

  endif

  n_points=c  

  deallocate(dat)

 end subroutine interpolate_1D_inv                 
 !******************************************************************       
 subroutine mask_for_2D_inv(targt,t2,datcomp,relevent_points)
 !******************************************************************
 ! finds the points bordering triangles where f=targt, sets them to 1
 ! all other points are set to zero
 ! treats negative values of t2%dat as invalid
 !need to tidy up a bit
 !****************************************************************** 
  use S95tables_type2
  use S95tables_type3
  implicit none
  !--------------------------------------input
  double precision, intent(in) :: targt
  integer, intent(in) :: datcomp
  type(table2), intent(in) :: t2
  !-------------------------------------output
  integer, intent(out) :: relevent_points(:,:)
  !-----------------------------------internal
  integer :: i,n,nn,j
  double precision :: Mj_middle,Mi_middle,Mj
  logical :: intable
  double precision :: c_x(3),c_y(3),c_f(3)
  integer :: c_xint(3),c_yint(3)
  !-------------------------------------------        

  if(ubound(t2%dat,dim=1).ne.ubound(relevent_points,dim=1))stop'wrong input to mask_for_2D_inv 1'
  if(ubound(t2%dat,dim=2).ne.ubound(relevent_points,dim=2))stop'wrong input to mask_for_2D_inv 2'

  if(lbound(relevent_points,dim=1).ne.1)stop'wrong input to mask_for_2D_inv 3'
  if(lbound(relevent_points,dim=2).ne.1)stop'wrong input to mask_for_2D_inv 4'

  relevent_points = 0

  do j=1,t2%nx2-1
    Mj_middle=t2%xmin2+t2%sep2*dble(j-1) + 0.5D0*t2%sep2
    do i=1,t2%nx1-1      
      Mi_middle=t2%xmin1+t2%sep1*dble(i-1) + 0.5D0*t2%sep1
      do n=1,2         
        Mj=Mj_middle+ (-1.0D0)**dble(n)*t2%sep2*0.25D0

        call findcorners_tabletype2(Mi_middle,Mj,t2,datcomp,intable,c_x,c_y,c_f)

        if(   ( targt.le. maxval( c_f ) ) .and. &
          &   ( targt.gt. minval( c_f ) ) .and. &
          &   ( minval( c_f ).ge. 0.0D0 ) )then

           do nn=1,3
              c_xint(nn)=nint( (c_x(nn)-t2%xmin1)/t2%sep1 ) + 1
              c_yint(nn)=nint( (c_y(nn)-t2%xmin2)/t2%sep2 ) + 1
              relevent_points(c_yint(nn),c_xint(nn))=1
           enddo
        endif

      enddo
    enddo
  enddo

 end subroutine mask_for_2D_inv 

end module interpolate
