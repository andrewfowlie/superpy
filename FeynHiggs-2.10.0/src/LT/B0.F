* B0.F
* the scalar two-point function
* this function is adapted from Ansgar Denner's bcanew.f
* to the conventions of LoopTools
* this file is part of FeynHiggs 
* last modified 30 Nov 11 th

#include "externals.h"
#include "types.h"


	ComplexType function B0(p, m1, m2)
	implicit none
	RealType p, m1, m2

#include "lt.h"

	ComplexType B0q
	external B0q

	B0 = B0q(p, m1, m2, mudim)
	end

************************************************************************

	ComplexType function B0q(p, m1, m2, q)
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
	    B0q = -log(m2/q) - fpv(1, x1, y1) - fpv(1, x2, y2)
	  else if( abs(x1) .lt. 10 .and. abs(x2) .lt. 10 ) then
	    B0q = 2 - log(p*oneMeps/q) +
     &        xlogx(-x1) + xlogx(-x2) - xlogx(y1) - xlogx(y2)
	  else if( abs(x1) .gt. .5D0 .and. abs(x2) .gt. .5D0 ) then
	    B0q = -log(m1/q) - fpv(1, y1, x1) - fpv(1, y2, x2)
	  else
	    print *, "B0(", p, ",", m1, ",", m2, ") not defined"
	    B0q = nan
	  endif
* zero momentum
	else if( abs(m1 - m2) .gt. minacc ) then
	  x2 = oneMeps*m1/(m1 - m2)
	  y2 = oneMeps*m2/(m2 - m1)
	  if( abs(y2) .gt. .5D0 ) then
	    B0q = -log(m2/q) - fpv(1, x2, y2)
	  else
	    B0q = -log(m1/q) - fpv(1, y2, x2)
	  endif
	else
	  B0q = -log(m2/q)
	endif
        B0q = B0q + delta
	end

