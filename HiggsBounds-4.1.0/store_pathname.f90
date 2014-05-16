!******************************************************************
module store_pathname
!******************************************************************
 implicit none

 integer,parameter:: pathname_length= 59
 character(len=pathname_length),parameter :: pathname= &
     &     "/home/andrew/WORKING/git-superpy/superpy/HiggsBounds-4.1.0" // &
     &     "/"

end module store_pathname
!******************************************************************
