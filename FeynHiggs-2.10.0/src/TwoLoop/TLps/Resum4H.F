* Resum4H.F
* resum the quartic Higgs coupling
* this file is part of FeynHiggs
* last modified 21 Nov 13 th

#include "externals.h"
#include "types.h"
#include "debug.h"


#define igs2 1
#define iht2 2
#define ilam 3

#define Nf 6

	subroutine Resum4H(dPhi22)
	implicit none
	RealType dPhi22

#include "TLps.h"

	RealType MX2, logMX2, logMT2, llog, lam12
	RealType x, gs2MT, ht2MT, y(3)
	RealType XtOS, XtMS
	external rge

	RealType CNf
	parameter (CNf = 2/3D0*Nf - 11)

	MX2 = sqrt(MSf2(1,3,3)*MSf2(2,3,3))
c	XtOS = Xf(3,3)
	XtOS = Re(Xtc)
	gs2MT = 4*pi*AlfasMT
	ht2MT = 4*pi*Alfat

* high scale
	logMX2 = log(MX2)
* low scale
	logMT2 = log(Mf2(tT,3))
	llog = logMX2 - logMT2

	XtMS = (1 + llog*(AlfasMT/pi - 3/(64*pi**2)*ht2MT))*XtOS

* initial values @ MX using analytically integrated 1L RGEs
	x = 1 - CNf/(16*pi**2)*gs2MT*llog
	y(igs2) = gs2MT/x
	y(iht2) = 2*(8 + CNf)*gs2MT*ht2MT/
     &    ((2*(8 + CNf)*gs2MT - 9*ht2MT)*x**(-8/CNf) + 9*ht2MT*x)
	x = XtMS**2/MX2
	x = x*(1 - 1/12D0*x)
	y(ilam) = 3/(8*pi**2)*y(iht2)**2*x

	if( debuglevel .gt. 4 ) then
	  DHIGGS "running from logMX2 =", logMX2 ENDL
	  DHIGGS "          to logMT2 =", logMT2 ENDL
	  DHIGGS "inputs:"		ENDL
	  DHIGGS "gs2(MT) =", gs2MT	ENDL
	  DHIGGS "ht2(MT) =", ht2MT	ENDL
	  DHIGGS "rough values @ MX:"	ENDL
	  DHIGGS "gs2(MX) =", y(igs2)	ENDL
	  DHIGGS "ht2(MX) =", y(iht2)	ENDL
	  DHIGGS "lam(MX) =", y(ilam)	ENDL
	endif

	call ODESolve(3, rge, y, logMX2, logMT2, 1D-5, -2D0)

	if( debuglevel .gt. 4 ) then
	  DHIGGS "rough values @ MT:"	ENDL
	  DHIGGS "[tmp] gs2(MT) =", y(igs2)	ENDL
	  DHIGGS "[tmp] ht2(MT) =", y(iht2)	ENDL
	  DHIGGS "lam(MT) =", y(ilam)	ENDL
	endif

	y(igs2) = gs2MT
	y(iht2) = ht2MT
	call ODESolve(3, rge, y, logMT2, logMX2, 1D-5, +2D0)

	if( debuglevel .gt. 4 ) then
	  DHIGGS "refined values @ MX:"	ENDL
	  DHIGGS "gs2(MX) =", y(igs2)	ENDL
	  DHIGGS "ht2(MX) =", y(iht2)	ENDL
	  DHIGGS "[tmp] lam(MX) =", y(ilam)	ENDL
	endif

	y(ilam) = 3/(8*pi**2)*y(iht2)**2*x
	call ODESolve(3, rge, y, logMX2, logMT2, 1D-5, -2D0)

* subtract leading logs
	x = XtOS**2/MX2
	lam12 = 3/(8*pi**2)*ht2MT**2*(x*(1 - 1/12D0*x) + llog)
	ht2MT = Mf2(tM2,3)/Mf2(tT,3)*ht2MT
	llog = logMX2 - log(Mf2(tM2,3))
	lam12 = lam12 + 3/(8*pi**2)*ht2MT**2*llog*(
     &    1/(32*pi**2)*ht2MT*(3*(llog + x) - 10) -
     &    1/(6*pi**2)*gs2MT*(3*(llog + x) - 2) )

	dPhi22 = vev**2/SB2*(y(ilam) - lam12)

	if( debuglevel .gt. 4 ) then
	  DHIGGS "final values @ MT:"	ENDL
	  DHIGGS "gs2(MT) =", y(igs2)	ENDL
	  DHIGGS "ht2(MT) =", y(iht2)	ENDL
	  DHIGGS "lam(MT) =", y(ilam)	ENDL
	  DHIGGS "dPhi22  =", dPhi22	ENDL
	endif
	end


************************************************************************

	subroutine rge(logQ2, y, dydx)
	implicit none
	RealType logQ2, y(*), dydx(*)

#include "const.h"

	RealType gs2, ht2, lam2, lam

	RealType c1L, c2L
	parameter (c1L = 1/(16*pi**2))
	parameter (c2L = c1L**2)
c	parameter (c2L = 0)

	gs2 = y(igs2)
	ht2 = y(iht2)
	lam = y(ilam)
	lam2 = lam**2

	dydx(igs2) = gs2**2*(
     &    c1L*(2/3D0*Nf - 11) +
     &    c2L*((38/3D0*Nf - 102)*gs2 - 2*ht2) )

	dydx(iht2) = ht2*(
     &    c1L*(9/2D0*ht2 - 8*gs2) +
     &    c2L*(6*ht2*(6*gs2 - 2*ht2 - lam) + 3/2D0*lam2 +
     &      (40/9D0*Nf - 404/3D0)*gs2**2) )

	dydx(ilam) = .5D0*(
     &    c1L*12*((lam - ht2)*(lam + ht2) + lam*ht2) +
     &    c2L*(ht2*(3*ht2*(20*ht2 - lam) + gs2*(80*lam - 64*ht2) -
     &      72*lam2) - 78*lam*lam2) )
	end

