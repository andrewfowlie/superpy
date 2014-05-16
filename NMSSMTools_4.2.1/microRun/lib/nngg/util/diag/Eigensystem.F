* Eigensystem.F
* diagonalization of a Hermitian n-by-n matrix using the Jacobi algorithm
* code adapted from the "Handbook" routines for complex A
* (Wilkinson, Reinsch: Handbook for Automatic Computation, p. 202)
* this file is part of FormCalc
* last modified 18 Jun 04 th

#include "config.h"


************************************************************************
** Eigensystem diagonalizes a Hermitian n-by-n matrix.
** Input: n, A = n-by-n matrix, Hermitian
** (only the upper triangle of A needs to be filled).
** Output: d = vector of eigenvalues, U = transformation matrix
** these fulfill diag(d) = U A U^+.

	subroutine Eigensystem(n, A, lda, d, U, ldu)
	implicit none
	integer n, lda, ldu
	double complex A(lda,*), U(ldu,*)
	double precision d(*)

	integer sweep, p, q, j
	double precision abssum, absA, threshold, delta, tau, t
	double complex s, x, y
	double precision Z(2,MAXMAT)

	double precision eps
	parameter (eps = 5D-16)

	do p = 1, n
	  d(p) = dble(A(p,p))
	  Z(1,p) = d(p)
	  Z(2,p) = 0
	enddo

	do p = 1, n
	  do q = 1, n
	    U(q,p) = 0
	  enddo
	  U(p,p) = 1
	enddo

	do sweep = 1, 50
	  abssum = 0
	  do p = 1, n - 1
	    do q = p + 1, n
	      abssum = abssum + abs(A(p,q))
	    enddo
	  enddo
	  if( abssum .lt. eps ) goto 1

	  threshold = 0
	  if( sweep .lt. 4 ) threshold = abssum/(5*n**2)

	  do p = 1, n - 1
	    do q = p + 1, n
	      absA = abs(A(p,q))
	      t = 100*absA
	      if( sweep .gt. 4 .and.
     &            abs(Z(1,p)) + t .eq. abs(Z(1,p)) .and.
     &            abs(Z(1,q)) + t .eq. abs(Z(1,q)) ) then
	        A(p,q) = 0
	      else
	        if( absA .gt. threshold ) then
	          call Jacobi(delta, s, tau,
     &              Z(1,p) - Z(1,q), A(p,q), absA)
	          Z(1,p) = Z(1,p) + delta
	          Z(2,p) = Z(2,p) + delta
	          Z(1,q) = Z(1,q) - delta
	          Z(2,q) = Z(2,q) - delta
	          A(p,q) = 0

	          do j = 1, p - 1
	            x = A(j,p)
	            y = A(j,q)
	            A(j,p) = x + dconjg(s)*y - tau*x
	            A(j,q) = y - s*x - tau*y
	          enddo

	          do j = p + 1, q - 1
	            x = A(p,j)
	            y = A(j,q)
	            A(p,j) = x + s*dconjg(y) - tau*x
	            A(j,q) = y - s*dconjg(x) - tau*y
	          enddo

	          do j = q + 1, n
	            x = A(p,j)
	            y = A(q,j)
	            A(p,j) = x + s*y - tau*x
	            A(q,j) = y - dconjg(s)*x - tau*y
	          enddo

	          do j = 1, n
	            x = U(p,j)
	            y = U(q,j)
	            U(p,j) = x + s*y - tau*x
	            U(q,j) = y - dconjg(s)*x - tau*y
	          enddo
	        endif
	      endif
	    enddo
	  enddo

	  do p = 1, n
	    d(p) = d(p) + Z(2,p)
	    Z(1,p) = d(p)
	    Z(2,p) = 0
	  enddo
	enddo

	print *, "Bad convergence in Eigensystem."

1	continue

* sort the eigenvalues

	do p = 1, n - 1
	  j = p
	  t = d(p)
	  do q = p + 1, n
#ifdef DESC
	    if( d(q) .gt. t ) then
#else
	    if( d(q) .lt. t ) then
#endif
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

