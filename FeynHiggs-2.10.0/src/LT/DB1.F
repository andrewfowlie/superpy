* DB1.F
* the derivative of B1 with respect to p^2
* this function is adapted from Ansgar Denner's bcanew.f
* to the conventions of LoopTools
* this file is part of FeynHiggs 
* last modified 30 Nov 11 th

#include "externals.h"
#include "types.h"


	ComplexType function DB1(p, m1, m2)
	implicit none
	RealType p, m1, m2

#include "lt.h"

	ComplexType fpv, yfpv
	external fpv, yfpv

	ComplexType x1, x2, y1, y2, r
	RealType minacc

	minacc = acc*(m1 + m2)
	if( abs(p) .gt. minacc ) then
* IR divergent case
	  if( m2 .eq. 0 .and. abs(p - m1) .lt. acc ) then
	    DB1 = .5D0*(3 + log(lambda/p))/p
	    return
	  endif
* general case
	  call roots(p, m1, m2, x1, x2, y1, y2, r)
	  if( abs(x1 - x2) .gt. acc*abs(x1 + x2) ) then
	    DB1 = (yfpv(2, x1, y1) - yfpv(2, x2, y2))/r
	  else if( abs(x1) .gt. 10 ) then
	    DB1 =
     &        Re((2/3D0 + (2 - 3*x1)*fpv(3, x1, y1))/x1**2)/p
	  else if( abs(y1) .gt. acc ) then
	    DB1 =
     &        Re(3/2D0 + (2 - 3*x1)*fpv(1, x1, y1))/p
	  else
	    print *, "DB1(", p, ",", m1, ",", m2, ") not defined"
	    DB1 = nan
	  endif
* zero momentum
	else if( abs(m1 - m2) .gt. minacc ) then
	  x2 = oneMeps*m1/(m1 - m2)
	  y2 = oneMeps*m2/(m2 - m1)
	  if( abs(x2) .lt. 10 ) then
	    DB1 = -(1/3D0 + yfpv(2, x2, y2))/(m1 - m2)
	  else
	    DB1 = -(1/3D0 + yfpv(3, x2, y2))/m1
	  endif
	else
	  DB1 = -1/12D0/m1
	endif
	end

