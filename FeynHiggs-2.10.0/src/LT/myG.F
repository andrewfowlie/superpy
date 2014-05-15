* myG.F
* this file is part of FeynHiggs
* last modified 30 Nov 11 th

#include "externals.h"
#include "types.h"


	RealType function myG(p, m1, m2, q)
	implicit none
	RealType p, m1, m2, q

	ComplexType A0q, B0q
	external A0q, B0q

	myG = (p - m1 - m2)*Re(B0q(p, m1, m2, q)) -
     &    Re(A0q(m1, q)) - Re(A0q(m2, q))
	end
