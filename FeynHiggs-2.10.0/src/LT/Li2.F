* Li2.F
* real dilogarithm
* this file is part of FeynHiggs
* last modified 26 Nov 13 th

#include "externals.h"
#include "types.h"


* ALGORITHM 490 COLLECTED ALGORITHMS FROM ACM
* ALGORITHM APPEARED IN [2]
* Real part of the dilogarithm function for a real argument
* [1] = L. Lewin, Dilogarithms and Associated Functions
* (Mac Donald, London, 1958).
* [2] = Comm. ACM, Vol. 18, No. 4N, p. 200

	RealType function Li2(x)
	implicit none
	RealType x

	RealType y, lx, ly, by, a, b, test
	integer n

	RealType pi, zeta2
	parameter (pi = 3.1415926535897932384626433832795029D0)
	parameter (zeta2 = pi**2/6D0)

* x0 = zero of dilog on the positive real axis
* lz = log(z)
	RealType x0, z, lz, c0, c1, c2, c3, c4
	parameter (x0 = 12.595170369845016128639896494864607D0)
	parameter (z = 1/(x0 - 1))
	parameter (lz = -2.4505886639885103433104596152498147D0)
	parameter (c0 = x0*z*.5D0)
	parameter (c1 = x0*z**2*(3*x0 - 2)/6D0)
	parameter (c2 = x0*z**3*((11*x0 - 15)*x0 + 6)/24D0)
	parameter (c3 = x0*z**4*
     &    (((50*x0 - 104)*x0 + 84)*x0 - 24)/120D0)
	parameter (c4 = x0*z**5*
     &    ((((274*x0 - 770)*x0 + 940)*x0 - 540)*x0 + 120)/720D0)

	if( x .gt. 12.6D0 ) goto 10
	if( x .ge. 12.59D0 ) goto 100
	if( x .ge. 2 ) goto 10
	if( x .gt. 1 ) goto 20
	if( x .eq. 1 ) goto 30
	if( x .gt. .5D0 ) goto 40
	if( x .gt. 1D-2 ) goto 50
	if( x .lt. -1 ) goto 60
	if( x .lt. -1D-2 ) goto 70

* dilog computed from [1] p. 244 Eq. (1)
	Li2 = x*(1 + x*(.25D0 + x*(1/9D0 + x*
     &    (625D-4 + x*(4D-2 + x*(1/36D0 + x*(1/49D0 +
     &    x/64D0)))))))
	return

* dilog computed from [1] p. 244 Eq. (6) and [2] Eq. (4)
10	y = 1/x
	by = -1 - y*(4 + y)
	Li2 = 2*zeta2 - .5D0*log(x)**2 +
     &    (y*(4 + 5.75D0*y) + 3*(1 + y)*(1 - y)*log(1 - y))/by
	if( Li2 + 4*y .eq. Li2 ) return
	goto 80

* dilog computed from [1] p. 244 Eq. (7) with x = 1/x + Eq. (6)
* and [2] Eq. (4)
20	y = 1 - 1/x
	lx = log(x)
	by = 1 + y*(4 + y)
	Li2 = zeta2 + lx*(.5D0*lx - log(x - 1)) +
     &    (y*(4 + 5.75D0*y) - 3*(1 + y)*lx/x)/by
	goto 80

* dilog computed from [1] p. 244 Eq. (2)
30	Li2 = zeta2
	return

* dilog computed from [1] p. 244 Eq. (7) and [2] Eq. (4)
40	y = 1 - x
	lx = log(x)
	by = -1 - y*(4 + y)
	Li2 = zeta2 - lx*log(y) +
     &    (y*(4 + 5.75D0*y) + 3*(1 + y)*lx*x)/by
	goto 80

* dilog computed from [2] Eq. (4)
50	y = x
	by = 1 + y*(4 + y)
	Li2 = (y*(4 + 5.75D0*y) + 3*(1 + y)*(1 - y)*log(1 - y))/by
	goto 80

* dilog computed from [1] p. 245 Eq. (12) with x = -x and [2] Eq. (4)
60	y = 1/(1 - x)
	lx = log(-x)
	ly = log(y)
	by = 1 + y*(4 + y)
	Li2 = -zeta2 + ly*(.5D0*ly + lx) +
     &    (y*(4 + 5.75D0*y) + 3*(1 + y)*(1 - y)*(lx + ly))/by
	if( Li2 + 4*y .eq. Li2 ) return
	goto 80

* dilog computed from [1] p. 244 Eq. (8) and [2] Eq. (4)
70	y = x/(x - 1)
	lx = log(1 - x)
	by = -1 - y*(4 + y)
	Li2 = -.5D0*lx**2 +
     &    (y*(4 + 5.75D0*y) - 3*(1 + y)*(1 - y)*lx)/by

80	b = 4*y**2/by
	do n = 1, 30
	  b = b*y
	  a = b/(n*(n + 1)*(n + 2))**2
	  test = Li2
	  Li2 = Li2 + a
	  if( Li2 .eq. test ) return
	enddo
	return

* dilog computed from Taylor series about zero of dilog, x0
100	y = x/x0 - 1
	a = 1 - y*(.5D0 - y*(1/3D0 - y*(.25D0 - y*(.2D0 - y/6D0))))
	b = c0 - y*(c1 - y*(c2 - y*(c3 - y*c4)))
	Li2 = y*(a*lz - y*b)
	end

