* SLHAClear.F
* clear (i.e. set to some bogus value) all data in the
* SUSY Les Houches Accord data structure (see hep-ph/0311123, 0801.0045)
* last modified 28 Jun 12 th

#include "SLHADefs.h"
#include "types.h"


	subroutine SLHAClear(slhadata)
	implicit none
	ComplexType slhadata(*)

	integer i

	do i = 1, OffsetDecays
	  slhadata(i) = invalid
	enddo

	do i = OffsetDecays + 1, nslhadata
	  slhadata(i) = ToComplex(invalid, invalid)
	enddo
	end

