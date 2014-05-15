* SLHAValid.F
* determine whether a block has all valid entries
* last modified 30 Nov 11 th

#include "SLHAIntern.h"


	integer function SLHAValid(slhablock, length)
	implicit none
	ComplexType slhablock(*)
	integer length

	integer i

	SLHAValid = 0
	do i = 1, length
	  if( Re(slhablock(i)) .eq. invalid ) return
	enddo
	SLHAValid = 1
	end

