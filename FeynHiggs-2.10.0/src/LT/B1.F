* B1.F
* the two-point tensor coefficient of p_mu
* this function is adapted from Ansgar Denner's bcanew.f
* to the conventions of LoopTools
* this file is part of FeynHiggs 
* last modified 30 Nov 11 th

#include "externals.h"
#include "types.h"


	ComplexType function B1(p, m1, m2)
	implicit none
	RealType p, m1, m2

#include "lt.h"

	ComplexType B1q
	external B1q

	B1 = B1q(p, m1, m2, mudim)
	end

************************************************************************

	ComplexType function B1q(p, m1, m2, q)
	implicit none
	RealType p, m1, m2, q

#include "lt.h"

	ComplexType fpv, xlogx
	external fpv, xlogx

	ComplexType x1, x2, y1, y2, r
	RealType minacc

	minacc = acc*(m1 + m2)
* general case
	if( abs(p) .gt. minacc ) then
	  call roots(p, m1, m2, x1, x2, y1, y2, r)
	  if( abs(y1) .gt. .5D0 .and. abs(y2) .gt. .5D0 ) then
	    B1q = .5D0*(log(m2/q) + fpv(2, x1, y1) + fpv(2, x2, y2))
	  else if( abs(x1) .lt. 10 .and. abs(x2) .lt. 10 ) then
	    B1q = -.5D0*(1 - log(p*oneMeps/q) +
     &        x1*xlogx(-x1) + x1 + x2*xlogx(-x2) + x2 -
     &        (1 + x1)*xlogx(y1) - (1 + x2)*xlogx(y2))
	  else if( abs(x1) .gt. .5D0 .and. abs(x2) .gt. .5D0 ) then
	    B1q = .5D0*(log(m1/q) + 1 +
     &        (1 + x1)*fpv(1, y1, x1) + (1 + x2)*fpv(1, y2, x2))
	  else
	    print *, "B1(", p, ",", m1, ",", m2, ") not defined"
	    B1q = nan
	  endif
* zero momentum
	else if( abs(m1 - m2) .gt. minacc ) then
	  x2 = oneMeps*m1/(m1 - m2)
	  y2 = oneMeps*m2/(m2 - m1)
	  if( abs(y2) .gt. .5D0 ) then
	    B1q = .5D0*(log(m2/q) + fpv(2, x2, y2))
	  else
	    B1q = .5D0*(log(m1/q) +
     &        (1 + x2)*fpv(1, y2, x2) + .5D0)
	  endif
	else
	  B1q = .5D0*log(m2/q)
	endif
	B1q = B1q - .5D0*delta
	end

