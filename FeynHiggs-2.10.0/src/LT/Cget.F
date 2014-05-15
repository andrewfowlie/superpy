* Cget.F
* the three-point tensor coefficients
* this file is part of FeynHiggs
* last modified 21 Sep 12 th

#include "externals.h"
#include "types.h"
#include "ltdefs.h"


	memindex function Cget(p1, p2, p1p2, m1, m2, m3)
	implicit none
	RealType p1, p2, p1p2, m1, m2, m3

#include "lt.h"
#include "ltcache.h"

	external Cijk

	RealType para(6)

	para(1) = p1
	para(2) = p2
	para(3) = p1p2
	para(4) = m1
	para(5) = m2
	para(6) = m3

	Cget = cacheindex(para, Cval(1,0), Cijk, 6, 13)
	end

************************************************************************

	subroutine Cput(res, p1, p2, p1p2, m1, m2, m3)
	implicit none
	ComplexType res(*)
	RealType p1, p2, p1p2, m1, m2, m3

#include "lt.h"
#include "ltcache.h"

	external Cijk

	RealType para(6)

	para(1) = p1
	para(2) = p2
	para(3) = p1p2
	para(4) = m1
	para(5) = m2
	para(6) = m3

	call cachecopy(res, para, Cval(1,0), Cijk, 6, Ncc)
	end

************************************************************************

	subroutine Cijk(C, P)
	implicit none
	ComplexType C(13)
	RealType P(6)

	ComplexType B0, B1, B11, C0
	external B0, B1, B11, C0

	RealType M11, M12, M22
	common /matrix2/ M11, M12, M22

	RealType det2, f1, f2
	ComplexType b023, b123, b1123
	ComplexType b013, b113, b1113
	ComplexType b012, b112, b1112
	ComplexType s1, s2

#define LC1() M11*s1 + M12*s2
#define LC2() M12*s1 + M22*s2

	M12 = .5D0*(P(2) - P(3) - P(1))
	det2 = 2*(P(1)*P(3) - M12*M12)
	M12 = M12/det2
	M11 = P(3)/det2
	M22 = P(1)/det2
	f1 = P(4) - P(5) + P(1)
	f2 = P(4) - P(6) + P(3)

	b023 = B0(P(2), P(6), P(5))
	b123 = B1(P(2), P(6), P(5))
	b1123 = B11(P(2), P(6), P(5))

	b013 = B0(P(3), P(4), P(6))
	b113 = B1(P(3), P(4), P(6))
	b1113 = B11(P(3), P(4), P(6))

	b012 = B0(P(1), P(4), P(5))
	b112 = B1(P(1), P(4), P(5))
	b1112 = B11(P(1), P(4), P(5))

	C(cc0) = C0(P(1), P(2), P(3), P(4), P(5), P(6))

	s1 = b013 - b023 - f1*C(cc0)
	s2 = b012 - b023 - f2*C(cc0)
	C(cc1) = LC1()
	C(cc2) = LC2()

	C(cc00) = .5D0*(P(4)*C(cc0) +
     &    .5D0*(f1*C(cc1) + f2*C(cc2) + b023)) + .25D0

	s1 = -f1*C(cc1) - b123 - 2*C(cc00)
	s2 = -f2*C(cc1) - b123 + b112
	C(cc11) = LC1()
	C(cc12) = LC2()

	b023 = b023 + b123
	s1 = b023 + b113 - f1*C(cc2)
	s2 = b023 - f2*C(cc2) - 2*C(cc00)
	C(cc12) = .5D0*(C(cc12) + LC1())
	C(cc22) = LC2()

	C(cc001) = (P(4)*C(cc1) + .5D0*(f1*C(cc11) +
     &    f2*C(cc12) + b123))/3D0 - 1/18D0
	C(cc002) = (P(4)*C(cc2) + .5D0*(f1*C(cc12) +
     &    f2*C(cc22) - b023))/3D0 - 1/18D0

	s1 = -b1123 - f1*C(cc11) - 4*C(cc001)
	s2 = -b1123 - f2*C(cc11) + b1112
	C(cc111) = LC1()
	C(cc112) = LC2()

	b1123 = b1123 + b023 + b123
	s1 = -b1123 - f1*C(cc22) + b1113
	s2 = -b1123 - f2*C(cc22) - 4*C(cc002)
	C(cc122) = LC1()
	C(cc222) = LC2()
	end

