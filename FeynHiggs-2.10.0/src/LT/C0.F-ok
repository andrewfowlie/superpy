* C0.F
* the scalar three-point function
* this function is adapted from Ansgar Denner's bcanew.f
* to the conventions of LoopTools
* this file is part of LoopTools
* last modified 1 Dec 11 th

#include "externals.h"
#include "types.h"


	ComplexType function C0(p1, p2, p1p2, m1, m2, m3)
	implicit none
	RealType p1, p2, p1p2, m1, m2, m3

#include "lt.h"

	ComplexType C0ir, C0reg
	external C0ir, C0reg

	if( m1 .eq. 0 .and.
     &      (abs(p1 - m2) + abs(p1p2 - m3)) .lt. acc ) then
	  C0 = C0ir(p2, p1, p1p2)
	else if( m2 .eq. 0 .and.
     &      (abs(p1 - m1) + abs(p2 - m3)) .lt. acc ) then
	  C0 = C0ir(p1p2, p1, p2)
	else if( m3 .eq. 0 .and.
     &      (abs(p2 - m2) + abs(p1p2 - m1)) .lt. acc ) then
	  C0 = C0ir(p1, p2, p1p2)
	else
	  C0 = C0reg(p1, p2, p1p2, m1, m2, m3)
	endif
	end

************************************************************************

	ComplexType function C0reg(p1, p2, p1p2, m1, m2, m3)
	implicit none
	RealType p1, p2, p1p2, m1, m2, m3

#include "lt.h"

	RealType q(5), m(5), mki, mkj, mij, qijk, ar
	ComplexType a, b, h, h0, h1, h2, h3, h4
	ComplexType y1, y2, y3, y4, x1, x2, x3, x4
	integer i, j, k

	ComplexType spence
	integer eta_n
	external spence, eta_n

	q(1) = p1
	q(2) = p2
	q(3) = p1p2
	q(4) = q(1)
	q(5) = q(2)

	m(1) = m1
	m(2) = m2
	m(3) = m3
	m(4) = m(1)
	m(5) = m(2)

	C0reg = 0

* all mom-squares != 0
	if( p1*p2*p1p2 .ne. 0 ) then
	  a = Sqrtc((p2 - p1 - p1p2)**2 - 4*p1*p1p2)
	  do i = 1, 3
	    j = i + 1
	    k = i + 2
	    mki = m(k) - m(i)
	    mkj = m(k) - m(j)
	    mij = m(i) - m(j)
	    qijk = q(i) - q(j) - q(k)
	    h2 = .5D0/(a*q(i))

	    h = q(i)*(qijk + mki + mkj) - mij*(q(j) - q(k))
	    y1 = h2*(h + a*(q(i) - mij))
	    y2 = h2*(h - a*(q(i) + mij))
	    b = Sqrtc((q(i) - mij)**2 - 4*q(i)*m(j))
	    y3 = h2*(h + a*b)
	    y4 = h2*(h - a*b)

	    h0 = q(i)*(q(j)*q(k) + qijk*m(k) + mki*mkj) -
     &        mij*(q(j)*mki - q(k)*mkj)
	    qijk = q(j) - q(k) - q(i)
	    h3 = h0 + q(j)*qijk*m(i) + q(k)*(q(k) - q(i) - q(j))*m(j)
	    if( abs(y3) .lt. abs(y4) ) then
	      y3 = h3/(a**2*q(i)*y4)
	    else
	      y4 = h3/(a**2*q(i)*y3)
	    endif
	    if( a*b .ne. 0 ) then
	      y3 = y3 + Ieps*abs(a*b*y3)/(a*b)
	      y4 = y4 - Ieps*abs(a*b*y4)/(a*b)
	    else
	      y3 = y3*onePeps
	      y4 = y4*oneMeps
	    endif

	    h1 = h2*(h - a*(q(i) - mij))
	    if( abs(y1) .lt. abs(h1) ) then
	      h3 = (q(i) - q(j) + mki)*(mkj*q(i) - mij*q(j)) +
     &          (q(i) - mij)*(q(j) - mkj)*q(k)
	      y1 = h3/(a**2*q(i)*h1)
	    endif
	    h1 = h2*(h + a*(q(i) + mij))
	    if( abs(y2) .lt. abs(h1) ) then
	      h3 = (q(i) - q(k) + mkj)*(mki*q(i) + mij*q(k)) +
     &          (q(i) + mij)*(q(k) - mki)*q(j)
	      y2 = h3/(a**2*q(i)*h1)
	    endif

	    C0reg = C0reg +
     &        spence(y2/y3, 0D0) + spence(y2/y4, 0D0) -
     &        spence(y1/y3, 0D0) - spence(y1/y4, 0D0)

	    if( Im(a) .ne. 0 ) then
	      h3 = abs(b)/b*Ieps
	      x1 = -.5D0*(q(i) - mij + b)/q(i) - h3
	      x2 = -.5D0*(q(i) - mij - b)/q(i) - h3
	      x3 = -.5D0*(-q(i) - mij + b)/q(i) - h3
	      x4 = -.5D0*(-q(i) - mij - b)/q(i) - h3
	      h3 = 1/y3
	      h4 = 1/y4
	      h = log(y1)*(eta_n(x1, x2) + eta_n(h3, h4) -
     &                     eta_n(x1, h3) - eta_n(x2, h4) ) -
     &            log(y2)*(eta_n(x3, x4) + eta_n(h3, h4) -
     &                     eta_n(x3, h3) - eta_n(x4, h4) ) +
     &            log(y3)*(eta_n(x1, h3) - eta_n(x3, h3)) +
     &            log(y4)*(eta_n(x2, h4) - eta_n(x4, h4))
	      if( Im(a) .gt. 0 .and. q(i) .lt. 0 )
     &          h = h - log(y1/y2)
	      C0reg = C0reg + c2ipi*h
	    endif
	  enddo
	  C0reg = C0reg/a
	  return
	endif

* one mom-square zero
	if( (p2*p1 + p1p2*p2 + p1*p1p2) .ne. 0 ) then
	  if( p1 .ne. 0 ) then
	    if( p2 .eq. 0 ) then
	      m(1) = m2
	      m(2) = m3
	      m(3) = m1
	      q(1) = p2
	      q(2) = p1p2
	      q(3) = p1
	    else
	      m(1) = m3
	      m(2) = m1
	      m(3) = m2
	      q(1) = p1p2
	      q(2) = p1
	      q(3) = p2
	    endif
	    m(4) = m(1)
	    m(5) = m(2)
	    q(4) = q(1)
	    q(5) = q(2)
	  endif
	  ar = q(2) - q(3)
	  do i = 2, 3
	    j = i + 1
	    k = i + 2
	    mki = m(k) - m(i)
	    mkj = m(k) - m(j)
	    mij = m(i) - m(j)
	    qijk = q(i) - q(j) - q(k)

	    if( i .eq. 2 ) then
	      y1 = 2*q(2)*(mki + ar)
	      y2 = 2*q(2)*mki
	    else
	      y1 = 2*q(3)*mkj
	      y2 = 2*q(3)*(mkj - ar)
	    endif
	    h = q(i)*(qijk + mki + mkj) - mij*(q(j) - q(k))
	    b = Sqrtc((q(i) - mij)**2 - 4*q(i)*m(j))
	    y3 = h + ar*b
	    y4 = h - ar*b

	    h0 = q(i)*(q(j)*q(k) + qijk*m(k) + mki*mkj) -
     &        mij*(q(j)*mki - q(k)*mkj)
	    h3 = h0 + q(j)*(q(j) - q(k) - q(i))*m(i) + 
     &                q(k)*(q(k) - q(i) - q(j))*m(j)
	    h3 = 4*h3*q(i)
	    if( abs(y3) .lt. abs(y4) ) then
	      y3 = h3/y4
	    else
	      y4 = h3/y3
	    endif
	    qijk = ar/q(i)
	    if( qijk .ne. 0 ) then
	      y3 = y3 + Ieps/qijk*abs(qijk*y3)
	      y4 = y4 - Ieps/qijk*abs(qijk*y4)
	    else
	      y3 = y3*onePeps
	      y4 = y4*oneMeps
	    endif

	    C0reg = C0reg +
     &        spence(y2/y3, 0D0) + spence(y2/y4, 0D0) -
     &        spence(y1/y3, 0D0) - spence(y1/y4, 0D0)
	  enddo
	  C0reg = C0reg/ar
	  return
	endif

* two mom-squares zero
	if( p1p2 .eq. 0 ) then
	  if( p2 .ne. 0 ) then
	    m(1) = m3
	    m(2) = m1
	    m(3) = m2
	    q(1) = p1p2
	    q(2) = p1
	    q(3) = p2
	  else
	    m(1) = m2
	    m(2) = m3
	    m(3) = m1
	    q(1) = p2
	    q(2) = p1p2
	    q(3) = p1
	  endif
	  m(4) = m(1)
	  m(5) = m(2)
	  q(4) = q(1)
	  q(5) = q(2)
	endif

	mki = m(2) - m(3)
	mkj = m(2) - m(1)
	mij = m(3) - m(1)

	if( m(2) .ne. m(3) ) then
	  y1 = -q(3) - mkj
	  y2 = -mkj
	  qijk = -mkj - q(3)*m(2)/mki
	  y3 = qijk - Ieps*sign(1D0, -q(3)/mki)*abs(qijk)
	  C0reg = C0reg + spence(y2/y3, 0D0) - spence(y1/y3, 0D0)
	endif

	b = Sqrtc((q(3) - mij)**2 - 4*q(3)*m(1))
	h = q(3)*(q(3) + mki + mkj)
	y1 = 2*q(3)*mkj
	y2 = 2*q(3)*(q(3) + mkj)
	y3 = h - q(3)*b
	y4 = h + q(3)*b
	h0 = 4*q(3)**2*(q(3)*m(2) + mki*mkj)
	if( abs(y3) .lt. abs(y4) ) then
	  y3 = h0/y4
	else
	  y4 = h0/y3
	endif
	y3 = y3 - Ieps*abs(y3)
	y4 = y4 + Ieps*abs(y4)

	C0reg = -(C0reg +
     &    spence(y2/y3, 0D0) + spence(y2/y4, 0D0) -
     &    spence(y1/y3, 0D0) - spence(y1/y4, 0D0))/q(3)
	end

************************************************************************

	ComplexType function C0ir(p2, p1, p1p2)
	implicit none
	RealType p2, p1, p1p2

#include "lt.h"

	ComplexType spence, ln
	external spence, ln

	RealType a, h1, h2, h3, ps
	ComplexType c

	if( abs(p1p2) .lt. acc .or. abs(p1) .lt. acc ) then
	  print *, "C0ir: mass singular case"
	  C0ir = nan
	  return
	endif

	if( abs(p2) .lt. acc ) then
	  if( abs(p1 - p1p2) .lt. acc ) then
	    C0ir = -.5D0*log(p1/lambda)/p1
	  else
	    C0ir = -.25D0*log(p1p2*p1/lambda**2)*
     &        log(p1/p1p2)/(p1 - p1p2)
	  endif
	  return
	endif

	ps = p2 - p1 - p1p2
	a = ps**2 - 4*p1*p1p2
	if( a .lt. 0 )
     &    print *, "C0ir: complex square root not implemented"
	a = sqrt(a)
	if( ps .le. 0 ) then
	  h1 = .5D0*(a - ps)
	else
	  h1 = -2*p1*p1p2/(a + ps)
	endif
	ps = p2 - p1 + p1p2
	if( ps .le. 0 ) then
	  h2 = .5D0*(a - ps)
	else
	  h2 = -2*p2*p1p2/(a + ps)
	endif
	ps = p2 + p1 - p1p2
	if( ps .le. 0 ) then
	  h3 = .5D0*(a - ps)
	else
	  h3 = -2*p1*p2/(a + ps)
	endif

	c = ln(-a/p2, -1D0)
	C0ir = (-pi**2/6D0 +
     &    spence(ToComplex(h2/a), -1D0) + spence(ToComplex(h3/a), -1D0) -
     &    .5D0*(ln(-h2/p2, -1D0)**2 + ln(-h3/p2, -1D0)**2) +
     &    .25D0*(ln(-p1/p2, -1D0)**2 + ln(-p1p2/p2, -1D0)**2) -
     &    c*(ln(-h1/p2, -1D0) - c) +
     &    ln(-lambda/p2, -1D0)*ln(h1/sqrt(p1*p1p2), 1D0))/a
	end

************************************************************************

	integer function eta_n(c1, c2)
	implicit none
	ComplexType c1, c2

	integer eta
	external eta

	eta_n = eta(c1, c2, 0D0, 0D0, 0D0)
	end

