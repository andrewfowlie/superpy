* SelectUZ.F
* select UHiggs or ZHiggs for the couplings and decay rates
* this file is part of FeynHiggs
* last modified 5 Jun 12 th

#include "externals.h"
#include "types.h"
#include "debug.h"


	subroutine FHSelectUZ(error, uzint_, uzext_, mfeff_)
	implicit none
	integer error, uzint_, uzext_, mfeff_

#include "FH.h"
#define __SUBROUTINE__ "FHSelectUZ"

	error = 0
	Check(uzint_, 0, 2, "uzint")
	Check(uzext_, 0, 2, "uzext")
	Check(mfeff_, 0, 2, "mfeff")
	if( error .ne. 0 ) return

	uzint = uzint_
	uzext = uzext_
	mfeff = mfeff_
	coup_valid = 0
	end

