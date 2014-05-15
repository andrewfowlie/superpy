* T134p.F
* the derivative of T134
* code taken from HDecay
* this file is part of FeynHiggs
* last modified 30 Nov 11 th

#include "externals.h"
#include "types.h"


	RealType function T134p(m1, m2, mu)
	implicit none
	RealType m1, m2, mu

	RealType x1, x2, l1, l2

	RealType pi, zeta2
	parameter (pi = 3.1415926535897932384626433832795029D0)
	parameter (zeta2 = pi**2/6D0)

	RealType Li2
	external Li2

	x1 = min(m1, m2)
	x2 = max(m1, m2)

	l1 = log(mu/x1)
	l1 = x1*(l1*(l1 + 3) + 3.5D0)
	if( abs(x2 - x1) .lt. 1D-12 ) then
	  T134p = 2*l1/mu
	else
	  l2 = log(mu/x2)
	  l2 = x2*(l2*(l2 + 3) + 3.5D0)
	  T134p = (l1 + l2 + (x2 - x1)*Li2(1 - x1/x2) -
     &      .5D0*x1*log(x1/x2)**2)/mu
	endif
	end

