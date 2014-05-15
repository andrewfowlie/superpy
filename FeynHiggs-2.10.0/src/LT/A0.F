* A0.F
* the scalar one-point function
* this function is adapted from Ansgar Denner's bcanew.f
* to the conventions of LoopTools
* this file is part of FeynHiggs 
* last modified 30 Nov 11 th

#include "externals.h"
#include "types.h"


	ComplexType function A0(m)
	implicit none
	RealType m

#include "lt.h"

	ComplexType A0q
	external A0q

	A0 = A0q(m, mudim)
	end

************************************************************************

	ComplexType function A0q(m, q)
	implicit none
	RealType m, q

#include "lt.h"

	if( m .eq. 0 ) then
	  A0q = 0
	else
	  A0q = m*(1 - log(m/q) + delta)
	endif
	end

