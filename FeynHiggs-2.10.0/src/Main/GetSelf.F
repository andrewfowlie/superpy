* GetSelf.F
* Subroutines for outside access to the renormalized Higgs self-energies
* this file is part of FeynHiggs
* last modified 30 Nov 11 th

#include "externals.h"
#include "types.h"
#include "debug.h"


	subroutine FHGetSelf(error, k2, key, sig, dkey, dsig)
	implicit none

#include "FH.h"
#define __SUBROUTINE__ "FHGetSelf"

	integer error
	RealType k2
	integer key, dkey
	ComplexType sig(semax), dsig(semax)

	integer i

	error = 0
	CheckSf(error)
	CheckTL(error)
	if( error .gt. 0 ) return

	if( higgs_valid .ne. valid ) call CalcCTs
	call CalcSelf(k2, key, dkey)

	do i = 1, semax
	  sig(i) = seR(i)
	  dsig(i) = dseR(i)
	enddo
	end
