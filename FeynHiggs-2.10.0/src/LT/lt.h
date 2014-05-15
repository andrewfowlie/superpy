* lt.h
* declarations internal to the LoopTools functions
* this file is part of FeynHiggs
* last modified 26 Nov 13 th


	RealType acc, pi, zeta2, eps
	ComplexType Ieps, onePeps, oneMeps, c2ipi, nan
	parameter (acc = 1D-12)
	parameter (pi = 3.1415926535897932384626433832795029D0)
	parameter (zeta2 = pi**2/6)
	parameter (eps = 1D-20)
	parameter (Ieps = (0D0,1D-20))
	parameter (onePeps = 1 + Ieps)
	parameter (oneMeps = 1 - Ieps)
	parameter (c2ipi = 2*pi*(0D0,1D0))
	parameter (nan = (1D123, 1D123))

	RealType mudim, delta, lambda
	common /cutoff/ mudim, delta, lambda

#ifndef ln
#define ln(x,s) log(x+(s)*Ieps)
#endif

