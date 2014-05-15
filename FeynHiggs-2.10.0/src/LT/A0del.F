* myA0del.F
* another version of the one-point function
* this file is part of FeynHiggs
* last modified 30 Nov 11 th

#include "externals.h"
#include "types.h"


	RealType function A0del(m)
	implicit none
	RealType m

#include "lt.h"

	RealType A0delq
	external A0delq

	A0del = A0delq(m, mudim)
	end

************************************************************************

	RealType function A0delq(m, q)
	implicit none
	RealType m, q

#include "lt.h"

	if( m .ne. 0 ) then
	  A0delq = .5D0*m*((1 - log(m/q))**2 + 1 + pi**2/6D0)
	else
	  A0delq = 0
	endif
	end

