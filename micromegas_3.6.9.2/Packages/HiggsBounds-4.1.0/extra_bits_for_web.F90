! This file is part of HiggsBounds
!  -KW
!******************************************************************
module extra_bits_for_web
! This files deals with the webversion of HiggsBounds
!******************************************************************

#ifdef NAGf90Fortran
 use F90_UNIX_ENV, only : iargc,getarg
 use F90_UNIX_IO, only : flush
#endif

 implicit none

 contains

 !************************************************************      
 subroutine getlongcommandline2web(d,gsq,partonicR,debugmode)
 !************************************************************
 ! reads components of d and dT from command line 
 ! note that debugmode for web version is different to other versions
 ! used for inputfile=F and whichinput='part' or 'hadr' only
 !************************************************************
  use usefulbits, only : np,Hneut,Hplus,Chineut,Chiplus,dataset,hadroncolliderextras,sqcouplratio 
  implicit none
  !--------------------------------------input
  type(dataset) :: d
  type(hadroncolliderextras) :: partonicR
  type(sqcouplratio) :: gsq
  logical :: debugmode
  !-----------------------------------internal
#ifndef NAGf90Fortran
  integer :: iargc
#endif
  double precision,allocatable :: CP_value_dble(:)
  integer :: nargs,i,j,k,n_iargc
  character(LEN=100) :: temp
  !-------------------------------------------   
  if(np(Hneut)>0)then
   allocate(CP_value_dble(np(Hneut)))   
  endif

  n_iargc=IARGC()

! the file "read_commandline_from_web.txt" has been generated using the same perl modules as the
! website, to try to ensure consistency 
#include "read_commandline_from_web.txt" 
  
  if( (k-1) .ne. nargs  )then
   write(*,*)'nargs=',nargs
   write(*,*)'k-1=',k-1
   call flush(6)
   stop "error in subroutine getlongcommandline2web (c)"
  endif    
   
  if(np(Hneut)>0)then
   d%CP_value=nint(CP_value_dble)
   deallocate(CP_value_dble)      
  endif

  contains

  subroutine nargs_errormsg
   if( n_iargc.ne. nargs  )then
    write(*,*)'nargs',nargs
    write(*,*)'IARGC()',IARGC()
    call flush(6)

    if(n_iargc.eq.5)then
     write(*,*) "At the moment, inputmethod=website"
     write(*,*) "Did you want inputmethod=datfile instead?"
     write(*,*) "(i.e. the command line version of HiggsBounds)" 
     write(*,*) "If so, please send the authors an email saying that"
     write(*,*) "they have left the website settings in accidently. Thanks!"
     write(*,*) "In the meantime, you can change inputmethod to"
     write(*,*) "inputmethod=datfile in the file HiggsBounds.F90"
     call flush(6)
     stop "Incorrect number of parameters given (getlongcommandline2web) A"
    else
     stop "Incorrect number of parameters given (getlongcommandline2web) B"
    endif

   endif
  end subroutine nargs_errormsg

 end subroutine getlongcommandline2web
 !********************************************************
 subroutine read_arg(x,k)
 !********************************************************  
  implicit none
  double precision :: x
  integer :: k
  character(LEN=100) :: temp
  
  temp=""  
  call GETARG(k,temp)
  read(temp,*)x   
  k=k+1
  
  end subroutine read_arg
 !********************************************************      
 subroutine weboutput(r)
 !********************************************************
 ! writes results to screen in a form suitable for the webpage
 ! used for inputfile=F only 
 !********************************************************
  use usefulbits, only : results,pr,listprocesses,vers
  use S95tables      
  implicit none
  !--------------------------------------input
  type(results) :: r
  !-----------------------------------internal  
  type(listprocesses) :: proc            
  integer :: ii
  character(LEN=200):: descrip 
  !-------------------------------------------      
             
  write(*,*)'********************************************************'
  write(*,*)
  if(abs(r%allowed95(1))==1)then
   write(*,*)'parameter point is UNEXCLUDED at 95 per cent C.L.'                     
  else
   write(*,*)'parameter point is EXCLUDED at 95 per cent C.L.'                   
  endif      

  write(*,*)'using the process with highest statistical sensitivity:'      
  
  proc=pr(r%chan(1))
  call outputproc(proc,6,descrip,1) 
  write(6,*)trim(descrip)
      
  write(*,*)'which has a theoretical rate vs. limit of'
  write(*,*)r%obsratio(1)  
      
  if(size(r%chan).gt.1)then      
      
   write(*,*)       
   write(*,*)'********************************************************'              
   write(*,*)'**************   Additional Information   **************'
   write(*,*)      
              
   do ii=2,size(r%chan)       
    select case(ii)
    case(2)
     write(*,*)'process with second highest statistical sensitivity:'
    case(3) 
     write(*,*)'process with third highest statistical sensitivity:'
    case default
     write(*,*)'another process:'
    end select    

    ! recall that 6 is standard output      
    proc=pr(r%chan(ii))
    call outputproc(proc,6,descrip,1) 
    write(6,*)trim(descrip)                            
    
    if(r%chan(ii).gt.0)then                   
      if(abs(r%allowed95(ii))==1)then
       write(*,*)'On its own, this process would have not excluded this'              
       write(*,*)'parameter point at 95 per cent C.L.'
      else
       write(*,*)'On its own, this process would have excluded this' 
       write(*,*)'parameter point at 95 per cent C.L.'
      endif
    
      write(*,*)'since it has a theoretical rate vs. limit of'
      write(*,*)r%obsratio(ii)                                  
      write(*,*)      
    endif     
   enddo
   
  endif            


  write(*,*)'using HiggsBounds version '//trim(adjustl(vers))//'.'
  write(*,*)'********************************************************'
            
            
  call webproclist

 end subroutine weboutput
 !********************************************************      
 subroutine webproclist
 ! writes the list of the processes considered by the webversion
 !********************************************************
  use usefulbits, only : pr,whichanalyses
  use S95tables, only : outputproc      
  implicit none
  integer :: x,y
  character(LEN=200):: descrip 

   do y=1,1
    write(*,*)
   enddo

   write(*,*)'List of processes considered in this calculation ('//trim(adjustl(whichanalyses))//'):'
   write(*,*)

   do x=1,ubound(pr,dim=1)
    if((pr(x)%findi.eq.1).and.(pr(x)%findj.eq.1))then
      call outputproc(pr(x),6,descrip,0)
      write(6,*)trim(descrip)       
    endif
   enddo   

 end subroutine webproclist
 !************************************************************
end module
!******************************************************************
