* DB11.F
* the derivative of B11 with respect to p^2
* this function is adapted from Ansgar Denner's bcanew.f
* to the conventions of LoopTools
* this file is part of FeynHiggs 
* last modified 30 Nov 11 th

#include "externals.h"
#include "types.h"


	ComplexType function DB11(p, m1, m2)
	implicit none
	RealType p, m1, m2

#include "lt.h"

	ComplexType fpv, yfpv
	external fpv, yfpv

	ComplexType x1, x2, y1, y2, r
	RealType minacc

	minacc = acc*(m1 + m2)
	DB11 = nan
* general case
	if( abs(p) .gt. minacc ) then
	  call roots(p, m1, m2, x1, x2, y1, y2, r)
	  if( abs(x1 - x2) .gt. acc*abs(x1 + x2) ) then
	    DB11 = (yfpv(3, x2, y2) - yfpv(3, x1, y1))/r
	  else if( abs(x1) .gt. 10 ) then
	    DB11 = Re((-3/4D0 + (4*x1 - 3)*fpv(4, x1, y1))/x1**2)/p
	  else if( abs(y1) .gt. acc ) then
	    DB11 = Re(-4/3D0 + (4*x1 - 3)*fpv(2, x1, y1))/p
	  endif
	else
* zero momentum
	  print *, "DB11(", p, ",", m1, ",", m2, ") not defined"
	endif
	end

