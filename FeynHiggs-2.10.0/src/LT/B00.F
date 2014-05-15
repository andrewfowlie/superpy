* B00.F
* the two-point tensor coefficient of g_{mu nu}
* this function is adapted from Ansgar Denner's bcanew.f
* to the conventions of LoopTools
* this file is part of FeynHiggs 
* last modified 30 Nov 11 th

#include "externals.h"
#include "types.h"


	ComplexType function B00(p, m1, m2)
	implicit none
	RealType p, m1, m2

#include "lt.h"

	ComplexType B00q
	external B00q

	B00 = B00q(p, m1, m2, mudim)
	end

************************************************************************

	ComplexType function B00q(p, m1, m2, q)
	implicit none
	RealType p, m1, m2, q

#include "lt.h"

	ComplexType A0q, B0q, B1q
	external A0q, B0q, B1q

	if( abs(p) .gt. acc*(abs(p) + m1 + m2) ) then
	  B00q = ((p + m1 - m2)*B1q(p, m1, m2, q) +
     &      2*m1*B0q(p, m1, m2, q) + A0q(m2, q) +
     &      m1 + m2 - p/3D0)/6D0
	else
	  B00q = (2*(m2*B0q(0D0, m1, m2, q) +
     &      A0q(m1, q)) + m1 + m2)/8D0
	endif
	end

