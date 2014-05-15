* DeltaRho_2L.F
* this file is part of FeynHiggs
* last modified 11 Apr 13 th

#include "externals.h"
#include "types.h"
#include "debug.h"


 	RealType function DeltaRho_2L_gluon()
        implicit none

#include "FH.h"

	ComplexType drho
	integer g

	ComplexType drf
	external drf

	drho = 0

	do g = 1, 3
	  drho = drho +
     &      USf2(1,1,3,g)*USf2(1,1,bBR,g)*
     &        drf(MSf2(1,3,g), MSf2(1,bBR,g)) +
     &      USf2(1,1,3,g)*USf2(2,1,bBR,g)*
     &        drf(MSf2(1,3,g), MSf2(2,bBR,g)) +
     &      USf2(2,1,3,g)*USf2(1,1,bBR,g)*
     &        drf(MSf2(2,3,g), MSf2(1,bBR,g)) +
     &      USf2(2,1,3,g)*USf2(2,1,bBR,g)*
     &        drf(MSf2(2,3,g), MSf2(2,bBR,g)) -
     &      USf2(1,2,3,g)*USf2(2,2,3,g)*
     &        drf(MSf2(1,3,g), MSf2(2,3,g)) -
     &      USf2(1,2,bBR,g)*USf2(2,2,bBR,g)*
     &        drf(MSf2(1,bBR,g), MSf2(2,bBR,g))
	enddo

	drho = 3*GF/(8*sqrt2*pi**2) * drho

	DeltaRho_2L_gluon = 2/(3*pi)*AlfasMT*Im(drho)

	if( debuglevel .lt. 3 ) return

	DCONST "DeltaRho(t,1L) =", 3*GF/(8*sqrt2*pi**2)*MT2	ENDL
	DCONST "DeltaRho(st|sb,1L,unused) =", Re(drho)		ENDL
	DCONST "DeltaRho(glu,1L+2L) =", DeltaRho_2L_gluon	ENDL

	if( debuglevel .lt. 4 ) return

	DCONST "DeltaRho: USf2(1,1,3,3) =", USf2(1,1,3,3)	ENDL
	DCONST "DeltaRho: USf2(1,2,3,3) =", USf2(1,2,3,3)	ENDL
	DCONST "DeltaRho: USf2(2,1,3,3) =", USf2(2,1,3,3)	ENDL
	DCONST "DeltaRho: USf2(2,2,3,3) =", USf2(2,2,3,3)	ENDL

	DCONST "DeltaRho: USf2(1,1,4,3) =", USf2(1,1,4,3)	ENDL
	DCONST "DeltaRho: USf2(1,2,4,3) =", USf2(1,2,4,3)	ENDL
	DCONST "DeltaRho: USf2(2,1,4,3) =", USf2(2,1,4,3)	ENDL
	DCONST "DeltaRho: USf2(2,2,4,3) =", USf2(2,2,4,3)	ENDL
	end


************************************************************************

	ComplexType function drf(x2, y2)
	implicit none
	RealType x2, y2

	ComplexType x, y, l, c
	RealType fy1, fs12

	ComplexType Dilog
	external Dilog

	ComplexType Ieps
	parameter (Ieps = (0D0, 1D-9))

	if( x2 .eq. y2 ) then
	  drf = 0
	  return
	endif

	x = (sqrt(x2) - Ieps)**2
	y = (sqrt(y2) - Ieps)**2
	l = x*log(x/y)
	c = l/(x - y)

	fy1 = Re(x + y - 2*y*c)
	fs12 = Re(x + y + c*(c*(x + y) - 2*l - 4*y) -
     &    2*(x - y)*Dilog(1 - x/y))

	drf = ToComplex(fy1, fs12)
	end

