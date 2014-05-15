* TakagiFactor.F
* computes the Takagi factorization of a complex symmetric matrix
* code adapted from the "Handbook" routines
* (Wilkinson, Reinsch: Handbook for Automatic Computation, p. 202)
* this file is part of FeynHiggs
* last modified 30 Nov 11 th

#include "externals.h"
#include "types.h"
#include "Util.h"


************************************************************************
** TakagiFactor factorizes a complex symmetric n-by-n matrix
** Input: n, A = n-by-n matrix, complex symmetric
** (only the upper triangle of A needs to be filled).
** Output: d = vector of diagonal values, U = transformation matrix
** these fulfill diag(d) = U^* A U^+ with U unitary.

	subroutine TakagiFactor(n, A,ldA, d, U,ldU, sort)
	implicit none
	integer n, ldA, ldU, sort
	ComplexType A(ldA,*), U(ldU,*)
	RealType d(*)

	integer p, q, j
	RealType red, off, thresh
	RealType sqp, sqq, t, invc
	ComplexType f, x, y
	ComplexType ev(2,MAXDIM)

	integer sweep
	common /nsweeps/ sweep

	if( n .gt. MAXDIM ) then
	  print *, "Dimension too large"
	  d(1) = -999
	  return
	endif

	do p = 1, n
	  ev(1,p) = 0
	  ev(2,p) = A(p,p)
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
	  if( sweep .lt. 4 ) thresh = off*red

	  do q = 2, n
	    do p = 1, q - 1
	      off = Sq(A(p,q))
	      sqp = Sq(ev(2,p))
	      sqq = Sq(ev(2,q))
	      if( sweep .gt. 4 .and.
     &            off .lt. SYM_EPS*(sqp + sqq) ) then
	        A(p,q) = 0
	      else if( off .gt. thresh ) then
	        t = .5D0*abs(sqp - sqq)
	        if( t .gt. 0 ) then
	          f = sign(1D0, sqp - sqq)*
     &              (ev(2,q)*Conjugate(A(p,q)) +
     &               Conjugate(ev(2,p))*A(p,q))
	        else
	          f = 1
	          if( sqp .ne. 0 ) f = sqrt(ev(2,q)/ev(2,p))
	        endif
	        t = t + sqrt(t**2 + Sq(f))
	        f = f/t

	        ev(1,p) = ev(1,p) + A(p,q)*Conjugate(f)
	        ev(2,p) = A(p,p) + ev(1,p)
	        ev(1,q) = ev(1,q) - A(p,q)*f
	        ev(2,q) = A(q,q) + ev(1,q)

	        t = Sq(f)
	        invc = sqrt(t + 1)
	        f = f/invc
	        t = t/(invc*(invc + 1))

	        do j = 1, p - 1
	          x = A(j,p)
	          y = A(j,q)
	          A(j,p) = x + (Conjugate(f)*y - t*x)
	          A(j,q) = y - (f*x + t*y)
	        enddo

	        do j = p + 1, q - 1
	          x = A(p,j)
	          y = A(j,q)
	          A(p,j) = x + (Conjugate(f)*y - t*x)
	          A(j,q) = y - (f*x + t*y)
	        enddo

	        do j = q + 1, n
	          x = A(p,j)
	          y = A(q,j)
	          A(p,j) = x + (Conjugate(f)*y - t*x)
	          A(q,j) = y - (f*x + t*y)
	        enddo

	        A(p,q) = 0

	        do j = 1, n
	          x = U(p,j)
	          y = U(q,j)
	          U(p,j) = x + (f*y - t*x)
	          U(q,j) = y - (Conjugate(f)*x + t*y)
	        enddo
	      endif
	    enddo
	  enddo

	  do p = 1, n
	    ev(1,p) = 0
	    A(p,p) = ev(2,p)
	  enddo
	enddo

	print *, "Bad convergence in TakagiFactor"

1	continue

* make the diagonal elements nonnegative

	do p = 1, n
	  d(p) = abs(A(p,p))
	  if( d(p) .gt. DBL_EPS .and. d(p) .ne. Re(A(p,p)) ) then
	    f = sqrt(A(p,p)/d(p))
	    do q = 1, n
	      U(p,q) = U(p,q)*f
	    enddo
	  endif
	enddo

	if( sort .eq. 0 ) return

* sort the eigenvalues

	do p = 1, n - 1
	  j = p
	  t = d(p)
	  do q = p + 1, n
	    if( sort*(t - d(q)) .gt. 0 ) then
	      j = q
	      t = d(q)
	    endif
	  enddo

	  if( j .ne. p ) then
	    d(j) = d(p)
	    d(p) = t
	    do q = 1, n
	      x = U(p,q)
	      U(p,q) = U(j,q)
	      U(j,q) = x
	    enddo
	  endif
	enddo
	end

