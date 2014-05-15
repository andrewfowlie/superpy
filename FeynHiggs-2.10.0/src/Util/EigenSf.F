* EigenSf.F
* diagonalization of the Hermitian 2-by-2 sfermion mass matrix
* this file is part of FeynHiggs
* last modified 5 Dec 13 th

#include "externals.h"
#include "types.h"
#include "debug.h"


	subroutine EigenSf(fail, ld, MSx, MSx2, U, M11, Mx, Xx, sort)
	implicit none
	integer fail, ld
	RealType MSx(ld,*), MSx2(ld,*)
	ComplexType U(ld,2,ld,*)
	RealType M11(*), Mx
	ComplexType Xx
	integer sort

#include "FH.h"

	RealType m1, m2, M12sq, delta, h, c
	ComplexType M12, s

* Note the conventions for better numerical stability:
* - M11(1) is the difference M11 - M22,
* - M11(2) is M22.

* MSx2(1) smaller sfermion mass
* MSx2(2) larger sfermion mass
* MSx2(3) = MSx2(2) - MSx2(1)
* MSx2(4) = MSx2(1) + MSx2(2)

	m2 = M11(2) + Mx**2
	m1 = m2 + M11(1)
	M12 = Mx*Conjugate(Xx)

	if( debuglevel .ge. 6 ) then
	  DPARA "EigenSf: Mx =", Mx	ENDL
	  DPARA "EigenSf: Xx =", Xx	ENDL
	  DPARA "EigenSf: M11 =", m1	ENDL
	  DPARA "EigenSf: M22 =", m2	ENDL
	  DPARA "EigenSf: M12 =", M12	ENDL
	endif

	M12sq = Sq(M12)
	h = .5D0*(M11(1) + sign(sqrt(M11(1)**2 + 4*M12sq), M11(1)))
	if( h .ne. 0 ) h = 1/h
	delta = h*M12sq
	c = 1/sqrt(delta*h + 1)
	s = h*c*M12

	MSx2(1,3) = M11(1) + 2*delta - 1D-15
	MSx2(1,4) = 2*m2 + M11(1)
	m1 = m1 + delta
	m2 = m2 - delta

	if( debuglevel .gt. 6 ) then
	  DPARA "EigenSf: m1 =", m1	ENDL
	  DPARA "EigenSf: m2 =", m2	ENDL
	endif

	if( sort .ne. 0 .and. m1 .gt. m2 ) then
	  MSx2(1,1) = m2
	  MSx2(1,2) = m1
	  U(1,1,1,1) = -Conjugate(s)
	  U(1,2,1,2) = s
	  U(1,1,1,2) = c
	  U(1,2,1,1) = c
	else
	  MSx2(1,1) = m1
	  MSx2(1,2) = m2
	  MSx2(1,3) = -MSx2(1,3)
	  U(1,1,1,1) = c
	  U(1,2,1,2) = c
	  U(1,1,1,2) = s
	  U(1,2,1,1) = -Conjugate(s)
	endif

	if( m1 .lt. 0 .or. m2 .lt. 0 ) fail = 1

	MSx2(1,1) = max(1D-100, MSx2(1,1))
	MSx(1,1) = sqrt(MSx2(1,1))
	MSx2(1,2) = max(2D-100, MSx2(1,2))
	MSx(1,2) = sqrt(MSx2(1,2))
	end

