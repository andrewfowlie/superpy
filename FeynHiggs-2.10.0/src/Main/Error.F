* Error.F
* print the error messages and return the error code
* this file is part of FeynHiggs
* last modified 27 Oct 11 th

#include "externals.h"
#include "debug.h"

#define MAXREPEAT 3
#define RESETTIME 10


	subroutine mj(error, lineno, routine, msg)
	implicit none
	integer error, lineno
	character*(*) routine, msg

* Note: 'time' is considered intrinsic by some compilers and
* external by others, so we give just the type:
	integer time

	integer hash, i, now

	integer n
	parameter (n = 32)
	integer tab(n,3)
	save tab
	data tab /n*0, n*0, n*1/

	if( lineno .gt. 0 ) error = lineno

* this is the DJB hash
	hash = 5381
	do i = 1, len(msg)
	  hash = hash + ishft(hash, 5) + ichar(msg(i:i))
	enddo
	do i = 1, len(routine)
	  hash = hash + ishft(hash, 5) + ichar(routine(i:i))
	enddo

	now = time()

	do i = 1, n
	  if( tab(i,1) .eq. 0 ) then
	    tab(i,1) = hash
	    tab(i,2) = now
	    goto 1
	  endif
	  if( hash .eq. tab(i,1) ) then
	    if( now - tab(i,2) .gt. RESETTIME ) tab(i,3) = 0
	    tab(i,2) = now
	    if( tab(i,3) .gt. MAXREPEAT ) return
	    tab(i,3) = tab(i,3) + 1
	    if( tab(i,3) .gt. MAXREPEAT ) then
100	      format(1A,$)
	      print 100, "Suppressing further "
	    endif
	    goto 1
	  endif
	enddo

1	continue

#if VT100
101	format(5A)
	if( lineno .gt. 0 ) then
	  print 101, BOLD//RED, routine, " error: ", msg, RESET
	else if( lineno .ne. -1 ) then
	  print 101, BOLD//BLUE, routine, " warning: ", msg, RESET
	else
	  print 101, routine, ": ", msg
	endif
#else
101	format(3A)
	if( lineno .gt. 0 ) then
	  print 101, routine, " error: ", msg
	else if( lineno .ne. -1 ) then
	  print 101, routine, " warning: ", msg
	else
	  print 101, routine, ": ", msg
	endif
#endif

	end

