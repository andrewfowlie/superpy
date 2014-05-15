* EDM.F
* computation of the electric dipole moments
* this file is part of FeynHiggs
* last modified 15 May 12 th

#include "externals.h"
#include "types.h"


************************************************************************
* hep-ph/9811202v2, Eq. (7)

	RealType function edmF(z)
	implicit none
	RealType z

#include "const.h"

	ComplexType delta, xp, xm

	RealType eps
	parameter (eps = 1D-10)

	ComplexType Dilog
	external Dilog

	if( abs(z - .25D0) .le. 1D-6 ) then
	   edmF = 2*(1 + log(z))
	else
	  delta = sqrt(1 - 4*(z - cI*eps))
	  xp = .5D0*(1 + delta)
	  xm = .5D0*(1 - delta)

	  edmF = 2 + log(z) +
     &      Re( z/delta*
     &         (log(z)*(log(-xp/xm) - log(-xm/xp))
     &          + 2*(Dilog(1/xp) - Dilog(1/xm))) )
	endif
	end


************************************************************************

	subroutine CalcEDMs(edmeTh, edmn, edmHg)
	implicit none
	RealType edmeTh, edmn, edmHg

#include "FH.h"

	integer t, ty, s, c, n, n2, Nc
	RealType edm(2,2:4)
	RealType aew, sum, sumC, pre, fac, ratio
	ComplexType gamma, h1a, h1b, h2

	RealType edmA, edmB, edmC, edmF
	external edmF

	integer g
	parameter (g = 1)

	RealType kappa(2:4), m(2:4), r(2:4)
	RealType Q(2:4)
	data Q /Qe, Qu, Qd/

#define T3(t) (iand(t, 1) - .5D0)

* hep-ph/9708456v3, Eq. (12)
	edmA(ratio) = .5D0/(ratio - 1)**2 *
     &    (3 - ratio + 2*log(ratio)/(1 - ratio))

* hep-ph/9708456v3, Eq. (13)
	edmB(ratio) = .5D0/(ratio - 1)**2 *
     &    (1 + ratio + 2*ratio*log(ratio)/(1 - ratio))

* hep-ph/9708456v3, Eq. (37)
	edmC(ratio) = 1/6D0/(ratio - 1)**2 *
     &    (-26 + 10*ratio + (2*ratio - 18)*log(ratio)/(1 - ratio))

	kappa(2) = Mf(2,g)/(sqrt2*MW*CB)
	kappa(3) = Mf(3,g)/(sqrt2*MW*SB)
	kappa(4) = Mf(4,g)/(sqrt2*MW*CB)

	m(2) = Mf(2,g)
	m(3) = 6D-3 * .35D0
	m(4) = 9.5D-3 * .35D0

	r(2) = TB
	r(3) = 1/TB
	r(4) = TB

	do t = 2, 4
	  edm(1,t) = 0
	  edm(2,t) = 0
	enddo

	aew = AlfaGF/(4*pi*SW2)

* Gluino contribution
* from hep-ph/9708456v3, Eqs. (19) and (36)

	do t = 3, 4
	  sum = 0
	  sumC = 0

	  ty = t + (bBR - 4)*isB(t)
	  do s = 1, 2
	    gamma = USfC(s,2,ty,g)*USf(s,1,ty,g)*M_3/abs(M_3)
	    pre = Im(gamma) * MGl/MSf2(s,ty,g)
	    ratio = MGl2/MSf2(s,ty,g)
	    sum = sum + pre*edmB(ratio)
	    sumC = sumC + pre*edmC(ratio)
	  enddo

	  edm(1,t) = edm(1,t) + Q(t)*2*AlfasMT/(3*pi)*sum
	  edm(2,t) = edm(2,t) + GSMT*AlfasMT/(4*pi)*sumC
	enddo

* Neutralino contribution
* from hep-ph/9708456v3, Eq. (22) and (38)

	do t = 2, 4
	  ty = t + (bBR - 4)*isB(t)
	  n2 = 3 + iand(t, 1)
	  sum = 0

	  do n = 1, 4
	    h1b = SW/CW*Conjugate(ZNeu(n,1))
	    h1a = -sqrt2*((Q(t) - T3(t))*h1b +
     &        T3(t)*Conjugate(ZNeu(n,2)))
	    h1b = sqrt2*Q(t)*h1b
	    h2 = kappa(t)*Conjugate(ZNeu(n,n2))

	    do s = 1, 2
	      gamma = (USf(s,1,ty,g)*h1a - USf(s,2,ty,g)*h2) *
     &          (USfC(s,2,ty,g)*h1b - USfC(s,1,ty,g)*h2)
	      pre = Im(gamma) * MNeu(n)/MSf2(s,ty,g)
	      ratio = MNeu2(n)/MSf2(s,ty,g)
	      sum = sum + pre*edmB(ratio)
	    enddo
	  enddo

	  edm(1,t) = edm(1,t) + Q(t)*aew*sum
	  edm(2,t) = edm(2,t) + GSMT*aew*sum
	enddo

* Chargino contribution (electron)
* from hep-ph/9708456v3, Eq. (29)

	sum = 0

	do c = 1, 2
	  gamma = kappa(2)*UChaC(c,2)*VChaC(c,1)
	  pre = Im(gamma) * MCha(c)/MSf2(1,1,g)
	  ratio = MCha2(c)/MSf2(1,1,g)
	  sum = sum + pre*Qe*edmA(ratio)
	enddo

	edm(1,2) = edm(1,2) + aew*sum

* Chargino contribution (up)
* from hep-ph/9708456v3, Eq. (27) and (39)

	sum = 0
	sumC = 0

	do s = 1, 2
	  do c = 1, 2
	    gamma = kappa(3)*VChaC(c,2)*USfC(s,1,bBR,g) *
     &        (UChaC(c,1)*USf(s,1,bBR,g) -
     &          kappa(4)*UChaC(c,2)*USf(s,2,bBR,g))
	    pre = Im(gamma) * MCha(c)/MSf2(s,bBR,g)
	    ratio = MCha2(c)/MSf2(s,bBR,g)
	    sum = sum + pre*(Qd*edmB(ratio) + (Qu - Qd)*edmA(ratio))
	    sumC = sumC + pre*edmB(ratio)
	  enddo
	enddo

	edm(1,3) = edm(1,3) + aew*sum
	edm(2,3) = edm(2,3) + GSMT*aew*sumC

* Chargino contribution (down)
* from hep-ph/9708456v3, Eq. (28)

	sum = 0
	sumC = 0

	do s = 1, 2
	  do c = 1, 2
	    gamma = kappa(4)*UChaC(c,2)*USfC(s,1,3,g) *
     &        (VChaC(c,1)*USf(s,1,3,g) -
     &          kappa(3)*VChaC(c,2)*USf(s,2,3,g))
	    pre = Im(gamma) * MCha(c)/MSf2(s,3,g)
	    ratio = MCha2(c)/MSf2(s,3,g)
	    sum = sum + pre*(Qu*edmB(ratio) + (Qd - Qu)*edmA(ratio))
	    sumC = sumC + pre*edmB(ratio)
	  enddo
	enddo

	edm(1,4) = edm(1,4) + aew*sum
	edm(2,4) = edm(2,4) + GSMT*aew*sumC

* two-loop contribution
* from hep-ph/0311314, Eq. (4.48) and hep-ph/9811202, Eq. (8)

	sum = 0
	sumC = 0

	do t = 2, 4
	  Nc = Ncolor(t)
	  ty = t + (bBR - 4)*isB(t)
	  fac = 2*Nc*Mf(t,3)*r(t)/(vev**2*SB*CB)
	  do s = 1, 2
	    pre = fac*edmF(MSf2(s,ty,3)/MA02)*
     &        Im(MUE * USf(s,1,ty,3) * USfC(s,2,ty,3))
	    sum = sum + Q(t)**2*pre
	    if( t .gt. 2 ) sumC = sumC + pre
	  enddo
	enddo

	do t = 2, 4
	  pre = r(t)*m(t)/MA02
	  edm(1,t) = edm(1,t) - AlfaGF/(32*pi**3)*Q(t)*pre*sum
	  edm(2,t) = edm(2,t) - GSMT*AlfasMT/(64*pi**3)*pre*sumC
	enddo

c FIXME: OVERALL SIGN OF CEDM ??

	edm(1,2) = GeV_cm*edm(1,2)
	edm(2,2) = GeV_cm*edm(2,2)
	edm(1,3) = GeV_cm/1.2D0*edm(1,3)
	edm(2,3) = GeV_cm/.9D0*edm(2,3)
	edm(1,4) = GeV_cm/1.2D0*edm(1,4)
	edm(2,4) = GeV_cm/.9D0*edm(2,4)

	edmeTh = edm(1,2)
	edmn = 2*.7D0*(edm(1,4) - .25D0*edm(1,3)) +
     &         2*.55D0*(edm(2,4) + .5D0*edm(2,3))
* extra factor 2 from Eq. (3.12) in hep-ph/0506106
	edmHg = .007D0*(edm(2,3) - edm(2,4)) + .01D0*edm(1,2)
	end

* current limits:
* |edmeTH| < 1.6 10^-27 e cm
* |edmHG|  < 2   10^-28 e cm
* |edmn|   < 6   10^-26 e cm

