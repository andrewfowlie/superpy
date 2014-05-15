* HiggsCorr.F
* Subroutines for the calculation of higher-order corrections
* to the MSSM Higgs sector
* this file is part of FeynHiggs
* V1.0 by Markus Frank 21/08/2002
* last modified 14 Sep 13 th

#include "externals.h"
#include "types.h"
#include "debug.h"

#define UNINITIALIZED -123D0


	subroutine FHHiggsCorr(error,
     &    MHiggs_, SAeff_, UHiggs_, ZHiggs_)
	implicit none

#include "FH.h"
#define __SUBROUTINE__ "FHHiggsCorr"

	integer error
	RealType MHiggs_(NHiggs)
	ComplexType SAeff_
	ComplexType UHiggs_(NNeutral,NNeutral)
	ComplexType ZHiggs_(NNeutral,NNeutral)

	integer i, j

	call CheckHiggsCorr(error)
	if( error .gt. 0 ) return

	do i = 1, NHiggs
	  MHiggs_(i) = MHiggs(i)
	enddo

	SAeff_ = SAeff

	do i = 1, NNeutral
	  do j = 1, NNeutral
	    UHiggs_(j,i) = UHiggs(j,i)
	    ZHiggs_(j,i) = ZHiggs(j,i)
	  enddo
	enddo
	end


************************************************************************
* main subroutine, calculates MHiggs, SAeff, and {U,Z}Higgs

	subroutine CheckHiggsCorr(error)
	implicit none
	integer error

#include "FH.h"

	integer i, j, k, l, p, uz, iter, minperm, hmix
	RealType m, d1, d2, test, mintest
	ComplexType p2, gw(NNeutral,NNeutral)
	ComplexType c, ZH(NNeutral,NNeutral,NNeutral)

	integer niter
	parameter (niter = 1)

	ComplexType GWFunc
	external GWFunc

	ComplexType CH(0:NNeutral,0:NNeutral,NNeutral)
	common /evmix/ CH

	integer perm(NNeutral,6), nperm(NNeutral)
	data perm /1,2,3, 2,1,3, 1,3,2, 2,3,1, 3,1,2, 3,2,1/
	data nperm /0, 2, 6/

	error = 0
	if( higgs_valid .eq. valid ) return
	CheckSf(error)
	coup_valid = 0

	entry CalcHiggsCorr(error)

	CheckTL(error)
	if( error .gt. 0 ) return

	call CalcCTs

	if( debuglevel .lt. 0 .and. higgs_valid .ne. valid ) then
	  call DumpHiggsFuncs
	  stop
	endif

	do uz = 0, 2
	  do i = 1, NNeutral
	    do j = 1, NNeutral
	      XHiggs(j,i,uz) = 0
	    enddo
	    XHiggs(i,i,uz) = 1
	  enddo
	enddo

	do i = 1, NHiggs
	  MHiggs2(i) = Msq(i)
	enddo

	hmix = higgsmix - ishft(higgsmix, -2)

	if( looplevel .ne. 0 ) then
	  call CalcUHiggs

	  do i = 1, NHiggs
	    if( debuglevel .ge. 3 )
     &        DHIGGS "searching Mh"//Digit(i) ENDL

	    p2 = ToComplex(MHiggs2(i), UNINITIALIZED)

	    if( i .eq. inputmass ) then
	      MHiggs2(i) = Msq(i)
	    else
	      call FindRoot(i, p2)
	      MHiggs2(i) = Re(p2)
	    endif

	    if( i .le. NNeutral ) then
	      k = ishft(hmix - i, -2)
	      do j = 1 + iand(k, i - 1),
     &               hmix + iand(k, i - hmix)
* The concept here is to distinguish the `stable' from the
* `unstable' zeros.  We do this by evaluating the GWFunc a
* little bit off the zero.  For the stable zero, the results
* are basically unchanged, but the unstable zero `rolls down'
* the hill and gives a completely different result.
	        c = 1.000003D0*p2
c	        c = p2 + 1D-1
	        do iter = 1, niter
	          gw(i,j) = GWFunc(j, c, ZH(1,i,j))
	        enddo
	      enddo
	    endif
	  enddo

	  if( debuglevel .ge. 3 ) then
	    DHIGGS "------------------------"	ENDL
	    DHIGGS "testing permutations:"	ENDL
	  endif

	  mintest = 1D300
	  minperm = 1
	  l = min(2, hmix)
	  do i = 1, nperm(hmix)
	    test = 0
	    do j = 1, hmix
	      p = perm(j,i)
	      d1 = 0
	      d2 = 0
	      do k = 1, hmix
	        d1 = d1 + abs(ZH(k,j,p) - CH(j,k,l))
	        d2 = d2 + abs(ZH(k,j,p) + CH(j,k,l))
	      enddo
	      test = test + abs(gw(j,p)) + min(d1, d2)
	    enddo
	    if( test .lt. mintest ) then
	      mintest = test
	      minperm = i
	    endif
	    if( debuglevel .ge. 3 )
     &        DHIGGS "test("//Digit(perm(1,i))//
     &                        Digit(perm(2,i))//
     &                        Digit(perm(3,i))//") =", test ENDL
	  enddo
	  if( debuglevel .ge. 3 )
     &      DHIGGS "minperm = "//Digit(perm(1,minperm))//
     &                           Digit(perm(2,minperm))//
     &                           Digit(perm(3,minperm)) ENDL
	  if( higgs_valid .ne. valid .and. mintest .gt. 25 )
     &      Warning("ZHiggs possibly unreliable")

	  do i = 1, NNeutral
	    p = perm(i,minperm)
	    do j = 1, NNeutral
	      ZHiggs(i,j) = ZH(j,i,p)
	    enddo
	  enddo
	endif

* sort the masses
	do i = 1, hmix - 1
	  p = i
	  m = MHiggs2(i)
	  do j = i + 1, hmix
	    if( MHiggs2(j) .lt. m ) then
	      m = MHiggs2(j)
	      p = j
	    endif
	  enddo

	  if( p .ne. i ) then
	    MHiggs2(p) = MHiggs2(i)
	    MHiggs2(i) = m
	    do j = 1, hmix
	      c = ZHiggs(i,j)
	      ZHiggs(i,j) = ZHiggs(p,j)
	      ZHiggs(p,j) = c
	    enddo
	  endif
	enddo

	do i = 1, NHiggs
	  MHiggs(i) = SignSqrt(MHiggs2(i))
	  if( MHiggs2(i) .lt. 0 ) then
* suppress the error message during FHUncertainties:
	    if( higgs_valid .eq. valid ) then
	      error = __LINE__
	    else
	      Error(error, "error computing Mh"//Digit(i))
	    endif
	  endif
	enddo

	if( hmix .le. 2 ) then
	  SAeff = SA*UHiggs(1,1) - CA*UHiggs(1,2)
	else
	  SAeff = SA
	endif

	higgs_valid = 0
	if( error .eq. 0 ) higgs_valid = valid

	if( debuglevel .lt. 2 ) return

	DHIGGS "------------------------"	ENDL
	DHIGGS "final values:"			ENDL
	DHIGGS "MHiggs(1) =", MHiggs(1)		ENDL
	DHIGGS "MHiggs(2) =", MHiggs(2)		ENDL
	DHIGGS "MHiggs(3) =", MHiggs(3)		ENDL
	DHIGGS "MHiggs(4) =", MHiggs(4)		ENDL
	DHIGGS "------------------------"	ENDL
	DHIGGS "UHiggs(1,1) =", UHiggs(1,1)	ENDL
	DHIGGS "UHiggs(1,2) =", UHiggs(1,2)	ENDL
	DHIGGS "UHiggs(1,3) =", UHiggs(1,3)	ENDL
	DHIGGS "UHiggs(2,1) =", UHiggs(2,1)	ENDL
	DHIGGS "UHiggs(2,2) =", UHiggs(2,2)	ENDL
	DHIGGS "UHiggs(2,3) =", UHiggs(2,3)	ENDL
	DHIGGS "UHiggs(3,1) =", UHiggs(3,1)	ENDL
	DHIGGS "UHiggs(3,2) =", UHiggs(3,2)	ENDL
	DHIGGS "UHiggs(3,3) =", UHiggs(3,3)	ENDL
	DHIGGS "------------------------"	ENDL
	DHIGGS "ZHiggs(1,1) =", ZHiggs(1,1)	ENDL
	DHIGGS "ZHiggs(1,2) =", ZHiggs(1,2)	ENDL
	DHIGGS "ZHiggs(1,3) =", ZHiggs(1,3)	ENDL
	DHIGGS "ZHiggs(2,1) =", ZHiggs(2,1)	ENDL
	DHIGGS "ZHiggs(2,2) =", ZHiggs(2,2)	ENDL
	DHIGGS "ZHiggs(2,3) =", ZHiggs(2,3)	ENDL
	DHIGGS "ZHiggs(3,1) =", ZHiggs(3,1)	ENDL
	DHIGGS "ZHiggs(3,2) =", ZHiggs(3,2)	ENDL
	DHIGGS "ZHiggs(3,3) =", ZHiggs(3,3)	ENDL
	DHIGGS "------------------------"	ENDL
	DHIGGS "SAeff =", SAeff			ENDL
	end


************************************************************************
* CalcUHiggs: Calculate UHiggs and initial values for MHiggs2 at k^2 = 0

	subroutine CalcUHiggs
	implicit none

#include "FH.h"

	integer key, lo
	ComplexType M(0:NNeutral,0:NNeutral)

* must cut off imaginary parts, otherwise UHiggs is not unitary
#define seU(i) Re(seR(i))

#define Need(i) -ibits(i-1-higgsmix,4,1)

	key = Key(h0h0) + Key(HHHH) + Key(HmHp) +
     &    iand(not(Key(inputmass)), Key(A0A0)) +
     &    iand(Need(2), Key(h0HH)) +
     &    iand(Need(3), Key(h0A0) + Key(HHA0)) +
     &    iand(Need(4), goldstones)

* calculate masses^2 & mixing angles at k^2 = 0
	call CalcSelf(-1D0, key, key)

	M(1,1) = Msq(h0h0) - seU(h0h0)
	M(2,2) = Msq(HHHH) - seU(HHHH)
	M(3,3) = Msq(A0A0) - seU(A0A0)
	M(1,2) = -seU(h0HH)
	M(1,3) = -seU(h0A0)
	M(2,3) = -seU(HHA0)

	M(0,0) = Msq(G0G0) - seU(G0G0)
	M(0,1) = -seU(h0G0)
	M(0,2) = -seU(HHG0)
	M(0,3) = -seU(A0G0)

	MHiggs2(HmHp) = Msq(HmHp) - seU(HmHp)

	if( debuglevel .ge. 4 ) then
	  DHIGGS "CalcUHiggs:"		ENDL
	  DHIGGS "M(1,1) =", M(1,1)	ENDL
	  DHIGGS "M(1,2) =", M(1,2)	ENDL
	  DHIGGS "M(2,2) =", M(2,2)	ENDL
	  DHIGGS "M(1,3) =", M(1,3)	ENDL
	  DHIGGS "M(2,3) =", M(2,3)	ENDL
	  DHIGGS "M(3,3) =", M(3,3)	ENDL
	  DHIGGS "M(0,0) =", M(0,0)	ENDL
	  DHIGGS "M(0,1) =", M(0,1)	ENDL
	  DHIGGS "M(0,2) =", M(0,2)	ENDL
	  DHIGGS "M(0,3) =", M(0,3)	ENDL
	endif

	lo = ieor(ishft(higgsmix, -2), 1)
	call HEigensystem(higgsmix, M(lo,lo),NNeutral+1, MHiggs2(lo),
     &    UHiggs(lo,lo),NNeutral+1, 1)

	if( debuglevel .lt. 3 ) return

	DHIGGS "========================"		ENDL
	DHIGGS "initial values:"			ENDL
	DHIGGS "MHiggs(0) =", SignSqrt(MHiggs2(0))	ENDL
	DHIGGS "MHiggs(1) =", SignSqrt(MHiggs2(1))	ENDL
	DHIGGS "MHiggs(2) =", SignSqrt(MHiggs2(2))	ENDL
	DHIGGS "MHiggs(3) =", SignSqrt(MHiggs2(3))	ENDL
	DHIGGS "MHiggs(4) =", SignSqrt(MHiggs2(4))	ENDL
	DHIGGS "------------------------"		ENDL
	end


************************************************************************
* EVFunc: the function whose roots are the exact poles of the
* Higgs propagators

	ComplexType function EVFunc(i, p2in)
	implicit none
	integer i
	ComplexType p2in

#include "FH.h"

	ComplexType CH(0:NNeutral,0:NNeutral,NNeutral)
	common /evmix/ CH

	integer key, lo, h
	ComplexType M(0:NNeutral,0:NNeutral), ev(0:NNeutral)
	ComplexType m2, p2

	RealType p2reim(2), p2re, p2im
	equivalence (p2, p2reim, p2re), (p2reim(2), p2im)

#define seEV(h) (seR(h) + cI*p2im*dseR(h))
*#define seEV(h) seR(h)

	p2 = p2in

	if( higgsmix - ishft(higgsmix, -2) .lt. max(2, i) ) then
	  key = iand(not(Key(inputmass)), Key(i))
	  call CalcSelf(p2re, key, key)
	  if( p2im .eq. UNINITIALIZED ) p2im = -Im(seR(i))

	  m2 = Msq(i) - seEV(i)
	else
	  key = cpeven +
     &      iand(Need(3), iand(not(Key(inputmass)), cpodd)) +
     &      iand(Need(4), goldstones)
	  call CalcSelf(p2re, key, key)
	  if( p2im .eq. UNINITIALIZED ) p2im = -Im(seR(i))

	  M(1,1) = Msq(h0h0) - seEV(h0h0)
	  M(2,2) = Msq(HHHH) - seEV(HHHH)
	  M(3,3) = Msq(A0A0) - seEV(A0A0)
	  M(1,2) = -seEV(h0HH)
	  M(1,3) = -seEV(h0A0)
	  M(2,3) = -seEV(HHA0)

	  M(0,0) = Msq(G0G0) - seEV(G0G0)
	  M(0,1) = -seEV(h0G0)
	  M(0,2) = -seEV(HHG0)
	  M(0,3) = -seEV(A0G0)

	  if( debuglevel .ge. 4 ) then
	    DHIGGS "EVFunc:"		ENDL
	    DHIGGS "M(1,1) =", M(1,1)	ENDL
	    DHIGGS "M(1,2) =", M(1,2)	ENDL
	    DHIGGS "M(2,2) =", M(2,2)	ENDL
	    DHIGGS "M(1,3) =", M(1,3)	ENDL
	    DHIGGS "M(2,3) =", M(2,3)	ENDL
	    DHIGGS "M(3,3) =", M(3,3)	ENDL
	    DHIGGS "M(0,0) =", M(0,0)	ENDL
	    DHIGGS "M(0,1) =", M(0,1)	ENDL
	    DHIGGS "M(0,2) =", M(0,2)	ENDL
	    DHIGGS "M(0,3) =", M(0,3)	ENDL
	  endif

	  lo = ieor(ishft(higgsmix, -2), 1)
	  call SEigensystem(higgsmix, M(lo,lo),NNeutral+1, ev(lo),
     &      CH(lo,lo,i),NNeutral+1, 1)
	  m2 = ev(i)

	  if( debuglevel .ge. 4 ) then
	    do h = lo, higgsmix
	      DHIGGS "ev("//Digit(h)//") =", sqrt(ev(h)) ENDL
	    enddo
	  endif
	endif

	p2in = ToComplex(Re(p2in), Im(m2))

c	p2 = sqrt(p2)
c	m2 = sqrt(m2)
	EVFunc = m2 - p2

	if( debuglevel .ge. 3 ) then
	  DHIGGS "         p =", p2 			ENDL
	  DHIGGS "  mh_ev("//Digit(i)//") =", m2	ENDL
	endif
	end


************************************************************************
* GWFunc: another function whose roots are the exact poles of the
* neutral Higgs propagators

	ComplexType function GWFunc(i, p2in, Z)
	implicit none

#include "FH.h"

	integer i
	ComplexType p2in, Z(NNeutral)

	integer j, k, ij, ik, jk
	integer needi, needj, needk, key
c	RealType mgamma
	ComplexType Dj, Dk, den, Zij, Zik, se_eff
	ComplexType dDj, dDk, dden, dZij, dZik, dse_eff
	ComplexType p2, m2

	RealType p2reim(2), p2re, p2im
	equivalence (p2, p2reim, p2re), (p2reim(2), p2im)

	integer se_ijk(NNeutral-1+NNeutral,NNeutral)
	data se_ijk /
     &    HHHH,A0A0, h0HH,h0A0,HHA0,
     &    h0h0,A0A0, h0HH,HHA0,h0A0,
     &    h0h0,HHHH, h0A0,HHA0,h0HH /

#define seGW(h) (seR(h) + cI*p2im*dseR(h))
*#define seGW(h) seR(h)
#define dseGW(h) dseR(h)

	p2 = p2in

	j = se_ijk(1,i)
	k = se_ijk(2,i)
	ij = se_ijk(3,i)
	ik = se_ijk(4,i)
	jk = se_ijk(5,i)

	needi = Need(i)
	needj = Need(j)
	needk = Need(k)
	key = iand( not(Key(inputmass)),
     &    Key(i) +
     &    iand(needj, Key(j) +
     &                iand(needi, Key(ij))) +
     &    iand(needk, Key(k) +
     &                iand(needi, Key(ik)) +
     &                iand(needj, Key(jk))) )
	call CalcSelf(p2re, key, key)
	if( p2im .eq. UNINITIALIZED ) p2im = -Im(seR(i))

	Dj = p2 - Msq(j) + seGW(j)
	Dk = p2 - Msq(k) + seGW(k)
	den = Dj*Dk - seGW(jk)**2
	if( den .eq. 0 ) then
* either j or k = inputmass and this boson does not mix
* (typically A0 and 2x2 mixing)
	  GWFunc = 1D100
	  return
	endif
	den = 1/den

	Zij = den*(seGW(ik)*seGW(jk) -
     &    seGW(ij)*Dk)
	Zik = den*(seGW(ij)*seGW(jk) -
     &    seGW(ik)*Dj)

	dDj = 1 + dseGW(j)
	dDk = 1 + dseGW(k)
	dden = dDj*Dk + Dj*dDk - 2*seGW(jk)*dseGW(jk)

	dZij = dseGW(ik)*seGW(jk) +
     &    seGW(ik)*dseGW(jk) -
     &    dseGW(ij)*Dk -
     &    seGW(ij)*dDk
	dZik = dseGW(ij)*seGW(jk) +
     &    seGW(ij)*dseGW(jk) -
     &    dseGW(ik)*Dj -
     &    seGW(ik)*dDj

	se_eff = seGW(i) +
     &    seGW(ij)*Zij +
     &    seGW(ik)*Zik
	dse_eff = dseGW(i) +
     &    dseGW(ij)*Zij +
     &    dseGW(ik)*Zik +
     &    den*(seGW(ij)*(dZij - dden*Zij) +
     &         seGW(ik)*(dZik - dden*Zik))

	if( debuglevel .ge. 4 ) then
	  DHIGGS "se_eff("//Digit(i)//") =", se_eff	ENDL
	  DHIGGS "dse_eff("//Digit(i)//") =", dse_eff	ENDL
	  DHIGGS "Z"//Digit(i)//Digit(j)//" =", Zij	ENDL
	  DHIGGS "Z"//Digit(i)//Digit(k)//" =", Zik	ENDL
	endif

	Z(i) = sqrt(1/(1 + dse_eff))
c version with jumps:
c	Z(i) = Sqrtc(1/(1 + Re(dse_eff)))
	Z(j) = Z(i)*Zij
	Z(k) = Z(i)*Zik

c	mgamma = Im(se_eff)/(1 + Re(dse_eff))
c	m2 = Msq(i) - Re(se_eff) - mgamma*Im(dse_eff)
c	m2 = Msq(i) - se_eff + cI*mgamma*dse_eff
	m2 = Msq(i) - se_eff

	p2in = m2

c	p2 = sqrt(p2)
c	m2 = sqrt(m2)
	GWFunc = m2 - p2

	if( debuglevel .ge. 3 ) then
c	  DHIGGS "  mh*width =", mgamma			ENDL
	  DHIGGS "         p =", p2 			ENDL
	  DHIGGS "  mh_gw("//Digit(i)//") =", m2	ENDL
	endif
	end


************************************************************************

	subroutine FindRoot(i, x)
	implicit none
	integer i
	ComplexType x

#include "FH.h"

	ComplexType lim(2,0:1), x1, x2, f1, f2, f
	RealType deriv_eps, s, m, w1, w2, s1, s2
	integer try, j

	ComplexType EVFunc
	external EVFunc

	RealType prec, inf, maxstep, widening_eps
	parameter (prec = 1D-3)
	parameter (inf = 2D0**1000)
	parameter (maxstep = 25)
	parameter (widening_eps = 25)
c	widening_eps = abs(Re(x))/2**17

	deriv_eps = abs(x)/2**27

* use Newton's method as long as the interval shrinks
	if( debuglevel .ge. 3 ) DHIGGS "FindRoot stage 1" ENDL
	lim(1,0) = inf
	lim(1,1) = -inf
	do try = 1, 20
	  x1 = x - deriv_eps
	  f1 = EVFunc(i, x1)
	  if( abs(f1) .lt. prec ) return
	  s = sign(1D0, Re(f1))
	  j = ishft(int(s) + 1, -1)

	  s = sign(1D0, Re(f1))
	  j = ishft(int(s) + 1, -1)
	  if( s*Re(lim(1,j) - x) .gt. widening_eps ) goto 2

	  x2 = x + deriv_eps
	  f2 = EVFunc(i, x2)

	  x = .5D0*(x1 + x2)
	  lim(1,j) = x
	  f = .5D0*(f1 + f2)
	  lim(2,j) = f

	  s = deriv_eps*Re(f)/Re(f2 - f1)
	  m = sqrt(abs(x))*maxstep
	  if( abs(s) .gt. m ) s = sign(m, s)
	  x = x - 2*s
	  if( abs(s) .lt. prec ) return
	enddo

2	if( abs(lim(1,1-j)) .eq. inf ) then
* find two points straddling the root
	  if( debuglevel .ge. 3 ) DHIGGS "FindRoot stage 2" ENDL
	  do try = 1, 20
	    s = Re(lim(2,j))
	    if( abs(s) .gt. maxstep ) s = sign(maxstep, s)
	    j = 1 - j
	    lim(1,j) = lim(1,1-j) + (2 + try)*s*sqrt(abs(lim(1,1-j)))
	    x = lim(1,j)
	    lim(2,j) = EVFunc(i, x)
	    if( Re(lim(2,j))*Re(lim(2,1-j)) .lt. 0 .and.
     &          Im(lim(2,j))*Im(lim(2,1-j)) .lt. 0 ) goto 3
	  enddo
	  return
	endif

* find the exact root using bisection
3	if( debuglevel .ge. 3 ) DHIGGS "FindRoot stage 3" ENDL
	do try = 1, 20
	  if( abs(Re(lim(1,1) - lim(1,0))) .lt. prec ) return
#if 1
	  w1 = 0
	  if( Re(lim(2,0))*Re(lim(2,1)) .lt. 0 ) then
	    s1 = abs(Re(lim(2,0))/(Re(lim(2,0)) - Re(lim(2,1))))
	    w1 = log(max(abs(Re(lim(2,0))), abs(Re(lim(2,1)))))
	  endif
	  w2 = 0
	  if( Im(lim(2,0))*Im(lim(2,1)) .lt. 0 ) then
	    s2 = abs(Im(lim(2,0))/(Im(lim(2,0)) - Im(lim(2,1))))
	    w2 = log(max(abs(Im(lim(2,0))), abs(Im(lim(2,1)))))
	  endif
	  s = .5D0
	  if( w1 + w2 .ne. 0 )
     &      s = min(max((w1*s1 + w2*s2)/(w1 + w2), .05D0), .95D0)
	  x = (1 - s)*lim(1,0) + s*lim(1,1)
#else
	  x = .5D0*lim(1,0) + .5D0*lim(1,1)
#endif
	  f = EVFunc(i, x)
	  if( abs(f) .lt. prec ) return

c	  s = sign(1D0, Re(f)*Re(lim(2,1)))
	  s = sign(1D0, Im(f)*Im(lim(2,1)))
	  j = ishft(int(sign(1D0, s)) + 1, -1)
	  lim(1,j) = x
	  lim(2,j) = f
	enddo
	end


************************************************************************

	subroutine DumpHiggsFuncs
	implicit none

#include "FH.h"

	integer pre_from, pre_to, pre_step
	integer pim_from, pim_to, pim_step
	integer pre, pim
	ComplexType p2, p2x, m2ev, m2gw, Zdummy(3)
	integer i, j, hfrom, hto
	character*20 env

	ComplexType EVFunc, GWFunc
	external EVFunc, GWFunc

	pre_from = 50
	pre_to = 200
	pre_step = 1

	pim_from = 0
	pim_to = 0
	pim_step = 5

	hfrom = 1
	hto = NHiggs
	if( debuglevel .ge. -4 ) then
	  hfrom = -debuglevel
	  hto = hfrom
	endif

	call getenv("FHPRERANGE", env)
	read(env, *, end=91, err=91) pre_from, pre_to, pre_step
91	continue

	call getenv("FHPIMRANGE", env)
	read(env, *, end=92, err=92) pim_from, pim_to, pim_step
92	continue

100	format("Writing {EV,GW}Func in range ", 3I8, " (re)"/
     &         "                             ", 3I8, " (im)")
	write(0,100) pre_from, pre_to, pre_step,
     &               pim_from, pim_to, pim_step

	do i = hfrom, hto
	  do pim = pim_from, pim_to, pim_step
	    do pre = pre_from, pre_to, pre_step
	      p2 = ToComplex(pre*abs(pre), pim*abs(pim))
* p2 gets overwritten by {EV,GW}Func
	      p2x = p2
	      do j = 1, 1
	        m2ev = EVFunc(i, p2)
	      enddo
	      if( i .eq. 4 ) then
	        m2gw = 0
	      else
	        m2gw = GWFunc(i, p2, Zdummy)
	      endif
101	      format(6G20.10)
	      write(70+i,101) p2x, m2ev, m2gw
	    enddo
102	    format(/)
	    write(70+i,102)
	  enddo
	enddo
	end

