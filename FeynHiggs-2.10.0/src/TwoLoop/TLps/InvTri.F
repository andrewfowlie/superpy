* InvTri.F
* the inverse of the triangle's area
* this file is part of FeynHiggs
* last modified 30 Nov 11 th

#include "externals.h"
#include "types.h"


	RealType function InvTri(a, b, c)
	implicit none
	RealType a, b, c

	InvTri = 1/(a**2 + b**2 + c**2 - 2*(a*b + a*c + b*c))
	end
