* GetSMPara.F
* get the derived SM parameters from the internal arrays
* this file is part of FeynHiggs
* last modified 30 Nov 11 th

#include "externals.h"
#include "types.h"
#include "debug.h"


	subroutine FHGetSMPara(error, CKM_)
	implicit none
	integer error
	ComplexType CKM_(3*3)

#include "FH.h"
#define __SUBROUTINE__ "FHGetSMPara"

	integer i

	error = 0
	CheckSMPara(error)

	do i = 1, 3*3
	  CKM_(i) = CKMin_flat(i)
	enddo
	end

