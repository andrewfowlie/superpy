* DB0.F
* the derivative of B0 with respect to p^2
* this function is adapted from Ansgar Denner's bcanew.f
* to the conventions of LoopTools
* this file is part of FeynHiggs 
* last modified 30 Nov 11 th

#include "externals.h"
#include "types.h"


	ComplexType function DB0(p, m1, m2)
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
	  if( m1*m2 .eq. 0 .and. abs(p - m1 - m2) .lt. acc ) then
	    DB0 = -(1 + .5D0*log(lambda/p))/p
	    return
	  endif
* general case
	  call roots(p, m1, m2, x1, x2, y1, y2, r)
	  if( abs(x1 - x2) .gt. acc*abs(x1 + x2) ) then
	    DB0 = (yfpv(1, x2, y2) - yfpv(1, x1, y1))/r
	  else if( abs(x1) .gt. 10 ) then
	    DB0 =
     &        Re(-(.5D0 + (1 - 2*x1)*fpv(2, x1, y1))/x1**2)/p
	  else if( abs(y1) .gt. acc .and. abs(x1) .gt. acc ) then
	    DB0 = -Re(2 + (1 - 2*x1)*fpv(0, x1, y1))/p
	  else
	    print *, "DB0(", p, ",", m1, ",", m2, ") not defined"
	    DB0 = nan
	  endif
* zero momentum
	else if( abs(m1 - m2) .gt. minacc ) then
	  x2 = oneMeps*m1/(m1 - m2)
	  y2 = oneMeps*m2/(m2 - m1)
	  if( abs(x2) .lt. 10 ) then
	    DB0 = (.5D0 + yfpv(1, x2, y2))/(m1 - m2)
	  else
	    DB0 = (.5D0 + yfpv(2, x2, y2))/m1
	  endif
	else
	  DB0 = 1/6D0/m1
	endif
	end

