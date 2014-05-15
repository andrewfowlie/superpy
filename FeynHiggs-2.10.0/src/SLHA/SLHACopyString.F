* SLHACopyString.F
* utility routine for SLHA{Put,Get}Info
* last modified 2 Oct 08 th

#include "SLHADefs.h"


	subroutine SLHACopyString(ns, s, nd, d)
	implicit none
	integer ns, nd
	character*1 s(ns), d(nd)

	integer m, i

	m = min(ns, nd)
	do i = 1, m
	  d(i) = s(i)
	enddo
	do i = m + 1, nd
	  d(i) = ' '
	enddo
	end

