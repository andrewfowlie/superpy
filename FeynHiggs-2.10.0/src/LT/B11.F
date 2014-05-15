* B11.F
* the two-point tensor coefficient of p_mu p_nu
* this function is adapted from Ansgar Denner's bcanew.f
* to the conventions of LoopTools
* this file is part of FeynHiggs 
* last modified 30 Nov 11 th

#include "externals.h"
#include "types.h"


	ComplexType function B11(p, m1, m2)
	implicit none
	RealType p, m1, m2

#include "lt.h"

	ComplexType B11q
	external B11q

	B11 = B11q(p, m1, m2, mudim)
	end

************************************************************************

	ComplexType function B11q(p, m1, m2, q)
	implicit none
	RealType p, m1, m2, q

#include "lt.h"

	ComplexType fpv, yfpv, xlogx
	external fpv, yfpv, xlogx

	ComplexType x1, x2, y1, y2, r
	RealType minacc

	minacc = acc*(m1 + m2)
* general case
	if( abs(p) .gt. minacc ) then
	  call roots(p, m1, m2, x1, x2, y1, y2, r)
	  if( abs(y1) .gt. .5D0 .and. abs(y2) .gt. .5D0 ) then
	    B11q = (-log(m2/q) - fpv(3, x1, y1) - fpv(3, x2, y2))/3D0
	  else if( abs(x1) .lt. 10 .and. abs(x2) .lt. 10 ) then
	    x1 = x1**2*(xlogx(-x1) + 1) + .5D0*x1 -
     &        (1 + x1*(1 + x1))*xlogx(y1)
	    x2 = x2**2*(xlogx(-x2) + 1) + .5D0*x2 -
     &        (1 + x2*(1 + x2))*xlogx(y2)
	    B11q = (2/3D0 - log(p*oneMeps/q) + x1 + x2)/3D0
	  else if( abs(x1) .gt. .5D0 .and. abs(x2) .gt. .5D0 ) then
	    x1 = (1 + x1*(1 + x1))*fpv(1, y1, x1) + .5D0*x1
	    x2 = (1 + x2*(1 + x2))*fpv(1, y2, x2) + .5D0*x2
	    B11q = (-4/3D0 - log(m1/q) - x1 - x2)/3D0
	  else
	    print *, "B11(", p, ",", m1, ",", m2, ") not defined"
	    B11q = nan
	  endif
* zero momentum
	else if( abs(m1 - m2) .gt. minacc ) then
	  x2 = oneMeps*m1/(m1 - m2)
	  y2 = oneMeps*m2/(m2 - m1)
	  if( abs(y2) .gt. .5D0 ) then
	    B11q = (-log(m2/q) - fpv(3, x2, y2))/3D0
	  else
	    B11q = (-log(m1/q) +
     &        (1 + x2*(1 + x2))*yfpv(0, x2, y2) +
     &        x2*(x2 + .5D0) + 1/3D0)/3D0
	  endif
	else
	  B11q = -log(m2/q)/3D0
	endif
	B11q = B11q + delta/3D0
	end

