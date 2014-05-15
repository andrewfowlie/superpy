* AddSelf.F
* Add user-specified values to the Higgs self-energies
* this file is part of FeynHiggs
* last modified 5 Nov 13 th

#include "externals.h"
#include "types.h"
#include "debug.h"


	subroutine FHAddSelf(error, sig, rotate)
	implicit none

#include "FH.h"
#define __SUBROUTINE__ "FHAddSelf"

	integer error, rotate
	ComplexType sig(semax)

	integer i

	error = 0
	CheckPara(error)
	if( error .gt. 0 ) return

	do i = 1, semax
	  seAdd(i) = sig(i)
	enddo

	if( rotate .ne. 0 ) then
	  seAdd(h0h0) = SA2*sig(h0h0) + CA2*sig(HHHH) - S2A*sig(h0HH)
	  seAdd(HHHH) = CA2*sig(h0h0) + SA2*sig(HHHH) + S2A*sig(h0HH)
	  seAdd(h0HH) = C2A*sig(h0HH) + .5D0*S2A*(sig(HHHH) - sig(h0h0))
	endif

	higgs_valid = 0
	coup_valid = 0

	if( debuglevel .lt. 3 ) return

	DSELF "User-defined additions to the self-energies:" ENDL
	if( rotate .ne. 0 ) then
	  DSELF "unrotated:"			ENDL
	  DSELF "sig(f1f1) +=", sig(h0h0)	ENDL
	  DSELF "sig(f2f2) +=", sig(HHHH)	ENDL
	  DSELF "sig(f1f2) +=", sig(h0HH)	ENDL
	  DSELF "rotated:"			ENDL
	endif
	DSELF "sig(h0h0) +=", seAdd(h0h0)	ENDL
	DSELF "sig(HHHH) +=", seAdd(HHHH)	ENDL
	DSELF "sig(h0HH) +=", seAdd(h0HH)	ENDL
	DSELF "sig(A0A0) +=", seAdd(A0A0)	ENDL
	DSELF "sig(HmHp) +=", seAdd(HmHp)	ENDL
	DSELF "sig(h0A0) +=", seAdd(h0A0)	ENDL
	DSELF "sig(HHA0) +=", seAdd(HHA0)	ENDL
	DSELF "sig(G0G0) +=", seAdd(G0G0)	ENDL
	DSELF "sig(h0G0) +=", seAdd(h0G0)	ENDL
	DSELF "sig(HHG0) +=", seAdd(HHG0)	ENDL
	DSELF "sig(A0G0) +=", seAdd(A0G0)	ENDL
	DSELF "sig(GmGp) +=", seAdd(GmGp)	ENDL
	DSELF "sig(HmGp) +=", seAdd(HmGp)	ENDL
	end

