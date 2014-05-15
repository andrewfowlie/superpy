* DB00.F
* the derivative of B00 with respect to p^2
* this function is adapted from Ansgar Denner's bcanew.f
* to the conventions of LoopTools
* this file is part of FeynHiggs 
* last modified 30 Nov 11 th

#include "externals.h"
#include "types.h"


	ComplexType function DB00(p, m1, m2)
	implicit none
	RealType p, m1, m2

#include "lt.h"

	ComplexType B1, DB0, DB1
	external B1, DB0, DB1

	DB00 = 1/6D0*(2*m1*DB0(p, m1, m2) +
     &    B1(p, m1, m2) +
     &    (p + m1 - m2)*DB1(p, m1, m2)) - 1/18D0
	end

