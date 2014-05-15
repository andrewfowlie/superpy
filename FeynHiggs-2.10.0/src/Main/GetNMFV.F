* GetNMFV.F
* get the derived NMFV parameters from the internal arrays
* this file is part of FeynHiggs
* last modified 30 Nov 11 th

#include "externals.h"
#include "types.h"
#include "debug.h"


	subroutine FHGetNMFV(error, MSS2_, Kf_)
	implicit none
	integer error
	ComplexType MSS2_(3*3*5), Kf_(3*3,2:4)

#include "FH.h"
#define __SUBROUTINE__ "FHGetNMFV"

	integer i, t

	error = 0
	CheckSf(error)
	if( error .gt. 0 ) return

	do i = 1, 3*3*5
	  MSS2_(i) = MSS2_flat(i)
	enddo

	do t = 2, 4
	  do i = 1, 3*3
	    Kf_(i,t) = Kfflat(i,t)
	  enddo
	enddo
	end

