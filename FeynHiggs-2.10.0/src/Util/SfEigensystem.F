* SfEigensystem.F
* diagonalization of a Hermitian n-by-n matrix using the Jacobi algorithm
* code adapted from the "Handbook" routines for complex A
* (Wilkinson, Reinsch: Handbook for Automatic Computation, p. 202)
* this file is part of FeynHiggs
* last modified 5 Jun 13 th

#include "externals.h"
#include "types.h"
#include "Util.h"


************************************************************************
** SfEigensystem diagonalizes a Hermitian n-by-n matrix.
** Input: n, A = n-by-n matrix, Hermitian
** (only the upper triangle of A needs to be filled).
** Output: d = vector of eigenvalues, U = transformation matrix
** these fulfill diag(d) = U A U^+ = U A U^-1 with U unitary.

	subroutine SfEigensystem(n, A,ldA, d, U,ldU, sort, stride)
	implicit none
	integer n, ldA, ldU, sort, stride
	ComplexType A(ldA,*), U(ldU,*)
	RealType d(*)

	integer p, q, pp, qq, j, o, r
	RealType red, off, thresh
	RealType t, delta, invc, s, App, Aqq
	ComplexType x, y, Apq
	RealType evp, evq, ev(2,MAXDIM)

	integer sweep
	common /nsweeps/ sweep

	if( n .gt. MAXDIM ) then
	  print *, "Dimension too large"
	  d(1) = -999
	  return
	endif

	do p = 1, n
	  ev(1,p) = 0
	  ev(2,p) = Re(A(p,p))
	  d(p) = ev(2,p)
	enddo

	do q = 2, n
	  do p = 1, q - 1
	    A(q,p) = Conjugate(A(p,q))
	  enddo
	enddo

	do p = 1, n
	  do q = 1, n
	    U(q,p) = 0
	  enddo
	  U(p,p) = 1
	enddo

	red = .04D0/n**4

	do sweep = 1, 50
	  off = 0
	  do q = 2, n
	    do p = 1, q - 1
	      off = off + Sq(A(p,q))
	    enddo
	  enddo
	  if( .not. off .gt. SYM_EPS ) goto 1

	  thresh = 0

	  do o = 0, n - 2
	    r = mod(stride - o, n) + (stride + o)/n*(n - 1)
	    do p = 1, n - r
	      q = p + r
	      Apq = A(p,q)
	      off = Sq(Apq)
	      if( sweep .gt. 4 .and. off .lt.
     &              SYM_EPS*(ev(2,p)**2 + ev(2,q)**2) ) then
	        A(p,q) = 0
	      else if( off .gt. thresh ) then
	        t = .5D0*(ev(2,p) - ev(2,q))
	        t = 1/(t + sign(sqrt(t**2 + off), t))

	        delta = t*off
	        invc = sqrt(delta*t + 1)
	        s = t/invc
	        t = delta/(invc + 1)

	        App = Re(A(p,p)) + s*(off - t*Re(A(p,p)))
	        Aqq = Re(A(q,q)) - s*(off + t*Re(A(q,q)))

	        if( o .ne. 0 .and. Re(s*Apq) .lt. 0 .and.
     &              abs(abs(Re(A(p,p)) - App) - abs(Re(A(p,p)) - Aqq) +
     &                  abs(Re(A(q,q)) - Aqq) - abs(Re(A(q,q)) - App))
     &                .lt. 1D-6 ) then
	          pp = q
	          qq = p
	          evq = d(q)
	          d(q) = d(p)
	          d(p) = evq
	        else
	          pp = p
	          qq = q
	        endif

	        evp = ev(1,p) + delta
	        evq = ev(1,q) - delta
	        ev(1,pp) = evp
	        ev(2,pp) = d(pp) + evp
	        ev(1,qq) = evq
	        ev(2,qq) = d(qq) + evq

	        do j = 1, n
	          x = A(j,p)
	          y = A(j,q)
	          A(j,pp) = x + s*(Conjugate(Apq)*y - t*x)
	          A(j,qq) = y - s*(Apq*x + t*y)
	        enddo

	        do j = 1, n
	          A(pp,j) = Conjugate(A(j,pp))
	          A(qq,j) = Conjugate(A(j,qq))
	        enddo

	        A(pp,pp) = App
	        A(qq,qq) = Aqq
	        A(pp,qq) = 0
	        A(qq,pp) = 0

	        do j = 1, n
	          x = U(p,j)
	          y = U(q,j)
	          U(pp,j) = x + s*(Apq*y - t*x)
	          U(qq,j) = y - s*(Conjugate(Apq)*x + t*y)
	        enddo
	      endif
	    enddo

	    if( sweep .lt. 4 ) thresh = off*red
	  enddo

	  do p = 1, n
	    ev(1,p) = 0
	    d(p) = ev(2,p)
	  enddo
	enddo

	print *, "Bad convergence in SfEigensystem"

1	if( sort .eq. 0 ) return

* sort the eigenvalues

	do p = 1, stride
	  q = p + stride
	  t = d(p)
	  if( sort*(t - d(q)) .gt. 0 ) then
	    d(p) = d(q)
	    d(q) = t
	    do j = 1, n
	      x = U(p,j)
	      U(p,j) = U(q,j)
	      U(q,j) = x
	    enddo
	  endif
	enddo
	end

