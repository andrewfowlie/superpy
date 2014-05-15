* SLHAExist.F
* determine whether a block has any valid entries
* last modified 30 Nov 11 th

#include "SLHAIntern.h"


	integer function SLHAExist(slhablock, length)
	implicit none
	ComplexType slhablock(*)
	integer length

	integer i

	SLHAExist = 0
	do i = 1, length
	  if( Re(slhablock(i)) .ne. invalid ) then
	    if( Im(slhablock(i)) .ne. 0 ) then
	      SLHAExist = 2
	      return
	    endif
	    SLHAExist = 1
	  endif
	enddo
	end

