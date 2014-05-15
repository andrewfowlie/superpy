* Phi.F
* this file is part of FeynHiggs
* last modified 1 Dec 11 th

#include "externals.h"
#include "types.h"


	RealType function Phi(x, y, z)
	implicit none
	RealType x, y, z

* from Davydychev and Tausk, Nucl. Phys. B397 (1993) 23

	RealType myPhi
	external myPhi

	if( x .le. z .and. y .le. z ) then
	  Phi = myPhi(x, y, z)
	else if( z .le. x .and. y .le. x ) then
	  Phi = z/x*myPhi(z, y, x)
	else
	  Phi = z/y*myPhi(z, x, y)
	endif
	end


************************************************************************
	
	RealType function myPhi(x, y, z)
	implicit none
	RealType x, y, z

	RealType u, v, upv, umv
	ComplexType clam, alt, cxp, cxm, ccphi

	ComplexType Dilog
	external Dilog

	RealType Pi
	parameter (Pi = 3.1415926535897932384626433832795029D0)

	u = sqrt(x/z)
	v = sqrt(y/z)
	upv = u + v
	umv = u - v
	clam = Sqrtc((1 - upv)*(1 + upv)*(1 - umv)*(1 + umv))

	cxp = .5D0*(1 + upv*umv - clam)
	alt = .5D0*(1 + upv*umv + clam)
	if( abs(cxp) .lt. abs(alt) ) cxp = u**2/alt
	cxm = .5D0*(1 - upv*umv - clam)
	alt = .5D0*(1 - upv*umv + clam)
	if( abs(cxm) .lt. abs(alt) ) cxm = v**2/alt

* phi function from Eq. (A4)

	ccphi = 2*(log(cxp)*log(cxm) - 2*log(u)*log(v) -
     &    (Dilog(cxp) + Dilog(cxm)) + Pi**2/6D0)/clam
	myPhi = Re(ccphi)
	end

