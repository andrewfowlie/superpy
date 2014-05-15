* auxB.F
* auxillary functions used by the two-point integrals
* these functions are adapted from Ansgar Denner's bcanew.f
* to the conventions of LoopTools
* this file is part of FeynHiggs 
* last modified 1 Dec 11 th

#include "externals.h"
#include "types.h"


	subroutine roots(p, m1, m2, x1, x2, y1, y2, r)
	implicit none
	RealType p, m1, m2
	ComplexType x1, x2, y1, y2, r

#include "lt.h"

	RealType qx, qy

	qx = p + m1 - m2
	qy = p - m1 + m2
	r = Sqrtc(p*(p - m1 - m2) - m1*qy - m2*qx)
	x1 = (qx + r)/2D0/p
	x2 = (qx - r)/2D0/p
	if( abs(x2) .gt. abs(x1) ) then
	  x1 = m1/(p*x2)
	else if( abs(x1) .gt. abs(x2) ) then
	  x2 = m1/(p*x1)
	endif
	x1 = x1 + abs(p*x1)/p*Ieps
	x2 = x2 - abs(p*x2)/p*Ieps

	y2 = (qy + r)/2D0/p
	y1 = (qy - r)/2D0/p
	if( abs(y2) .gt. abs(y1) ) then
	  y1 = m2/(p*y2)
	else if( abs(y1) .gt. abs(y2) ) then
	  y2 = m2/(p*y1)
	endif
	y1 = y1 - abs(p*y1)/p*Ieps
	y2 = y2 + abs(p*y2)/p*Ieps
	end

************************************************************************

	ComplexType function fpv(n, x, y)
	implicit none
	integer n
	ComplexType x, y

#include "lt.h"

	integer m
	ComplexType xm

	if( abs(x) .lt. 5 ) then
	  if( n .eq. 0 ) then
	    fpv = -log(-y/x)
	  else if( abs(x) .lt. acc ) then
	    fpv = -1D0/n
	  else
	    xm = -log(-y/x)
	    do m = 1, n
	      xm = x*xm - 1D0/m
	    enddo
	    fpv = xm
	  endif
	else
	  fpv = 0
	  xm = 1
	  do m = 1, 50
	    xm = xm/x
	    fpv = fpv + xm/(m + n)
	    if( abs(xm/fpv) .lt. acc**2 ) return
	  enddo
	endif
	end

************************************************************************

	ComplexType function yfpv(n, x, y)
	implicit none
	integer n
	ComplexType x, y

	ComplexType fpv
	external fpv

	if( abs(y) .eq. 0 ) then
	  yfpv = 0
	else
	  yfpv = y*fpv(n, x, y)
	endif
	end

************************************************************************

	ComplexType function xlogx(x)
	implicit none
	ComplexType x

	if( abs(x) .eq. 0 ) then
	  xlogx = 0
	else
	  xlogx = x*log(x)
	endif
	end

