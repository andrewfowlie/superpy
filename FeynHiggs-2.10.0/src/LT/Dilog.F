* Dilog.F
* complex dilogarithm
* this file is part of FeynHiggs
* last modified 2 Aug 12 th

#include "externals.h"
#include "types.h"


	ComplexType function Dilog(z)
	implicit none
	ComplexType z

	ComplexType Dilogsum
	external Dilogsum

	RealType absz, abs1z
	ComplexType t, mlogz

	RealType pi, zeta2
	parameter (pi = 3.1415926535897932384626433832795029D0)
	parameter (zeta2 = pi*pi/6D0)

	absz = abs(z)
	if( absz .lt. 1D-20 ) then
	  Dilog = -log(1 - z)
	  return
	endif

	abs1z = abs(1 - z)
	if( abs1z .lt. 1D-14 ) then
	  Dilog = zeta2
	  return
	endif

	if( Re(z) .gt. .5D0 ) then
	  mlogz = -log(z)
	  t = zeta2 + mlogz*log(1 - z)
	  if( abs1z .gt. 1 ) then
	    Dilog = Dilogsum(log(1 - 1/z)) + zeta2 + .5D0*log(z - 1)**2 + t
	  else
	    Dilog = -Dilogsum(mlogz) + t
	  endif
	else
	  if( absz .gt. 1 ) then
	    Dilog = -Dilogsum(-log(1 - 1/z)) - zeta2 - .5D0*log(-z)**2
	  else
	    Dilog = Dilogsum(-log(1 - z))
	  endif
	endif
	end


************************************************************************

	ComplexType function Dilogsum(w)
	implicit none
	ComplexType w

	ComplexType u, t
	integer k

	RealType b2, b4, b6, b8, b10, b12, b14
	RealType b16, b18, b20, b22
	parameter (b2 = 1/6D0)
	parameter (b4 = -1/30D0)
	parameter (b6 = 1/42D0)
	parameter (b8 = -1/30D0)
	parameter (b10 = 5/66D0)
	parameter (b12 = -691/2730D0)
	parameter (b14 = 7/6D0)
	parameter (b16 = -3617/510D0)
	parameter (b18 = 43867/798D0)
	parameter (b20 = -174611/330D0)
	parameter (b22 = 854513/138D0)

	RealType bernoulliB(11)
	data bernoulliB /b2, b4, b6, b8, b10, b12, b14,
     &    b16, b18, b20, b22/

	Dilogsum = w*(1 - .25D0*w)
	if( abs(w) .lt. 1D-10 ) return

	u = w
	do k = 1, 11
	  u = u*w**2/Re(2*k*(2*k + 1))
	  t = u*bernoulliB(k)
	  Dilogsum = Dilogsum + t
	  if( abs(t) .lt. 1D-16*abs(Dilogsum) ) return
	enddo
	end

