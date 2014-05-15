* TF.F
* this file is part of FeynHiggs
* last modified 1 Dec 11 th

#include "externals.h"
#include "types.h"


	ComplexType function TF(m1, m2, m3)
	implicit none
	RealType m1, m2, m3

	ComplexType Dilog
	external Dilog

	RealType m13, m23
	ComplexType lambda, a1, a2

	RealType pi, mudim
	parameter (pi = 3.1415926535897932384626433832795029D0)
	parameter (mudim = 1)

	m23 = m2/m3

	if( m1 .eq. 0 ) then
	  TF = (m3 - m2)*
     &      (Dilog(ToComplex(1 - m23)) - .5D0*log(m3/mudim)*log(m23))
c	  TF = (m3 - m2)*Dilog(ToComplex(1 - m23)) -
c     &      .5D0*log(m3/mudim)*log(m23)
	else
	  m13 = m1/m3
	  lambda = Sqrtc((m13 - m23)**2 - 2*(m13 + m23) + 1)
	  a1 = .5D0*(1 + m13 - m23 - lambda)
	  a2 = .5D0*(1 - m13 + m23 - lambda)
	  TF = m3*lambda*(log(a1)*log(a2) - Dilog(a1) - Dilog(a2) -
     &      .5D0*log(m13)*log(m23) + pi**2/6D0)
	endif
	end

