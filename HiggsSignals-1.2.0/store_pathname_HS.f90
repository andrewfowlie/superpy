!******************************************************************
module store_pathname_HS
!******************************************************************
 implicit none

 integer,parameter:: pathname_length= 60
 character(len=pathname_length),parameter :: pathname_HS= &
     &     "/home/andrew/WORKING/git-superpy/superpy/HiggsSignals-1.2.0" // &
     &     "/"

end module store_pathname_HS
!******************************************************************
