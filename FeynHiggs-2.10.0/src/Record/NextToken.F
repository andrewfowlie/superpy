* NextToken.F
* tokenize a string
* this file is part of FeynHiggs
* last modified 2 May 07 th

#include "externals.h"


	integer function NextToken(s, pos)
	implicit none
	character*(*) s
	integer pos

	integer c

	NextToken = 0

	do pos = pos, len(s)
	  c = ichar(s(pos:pos))
	  if( c .eq. 35 ) return
	  if( c .eq. 32 .or. c .eq. 9 .or. c .eq. 61 ) then
	    if( NextToken .ne. 0 ) return
	  else
	    if( NextToken .eq. 0 ) NextToken = pos
	  endif
	enddo
	end

