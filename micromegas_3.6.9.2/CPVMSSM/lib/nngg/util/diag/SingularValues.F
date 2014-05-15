* SingularValues.F
* singular-value decomposition of an m-by-n matrix
* this file is part of FormCalc
* last modified 18 Jun 04 th

#include "config.h"


************************************************************************
** SingularValues performs a singular-value decomposition.
** Input: m, n, A = m-by-n matrix (m > n).
** Output: d = vector of singular values,
** U = n-by-m left transformation matrix,
** V = n-by-n right transformation matrix,
** these fulfill diag(d) = U^* A V^+.

	subroutine SingularValues(m, n, A, lda, d, U, ldu, V, ldv)
	implicit none
	integer m, n, lda, ldu, ldv
	double complex A(lda,*)
	double precision d(*)
	double complex U(ldu,*), V(ldv,*)

	integer i, j, k, l, r
	double complex s
	double precision x

	double complex AdagA(MAXMAT,MAXMAT)

	double precision eps
	parameter (eps = 5D-16)

* The singular values are the roots of the eigenvalues of A^+ A

	do k = 1, n
	  do i = 1, k
	    s = 0
	    do j = 1, m
	      s = s + dconjg(A(j,i))*A(j,k)
	    enddo
	    AdagA(i,k) = s
	  enddo
	enddo
	call Eigensystem(n, AdagA, MAXMAT, d, V, ldv)

#ifdef DESC
	do k = 1, n
#else
	do k = n, 1, -1
#endif
	  if( d(k) .lt. eps ) goto 10
	  d(k) = sqrt(d(k))
	  x = 1/d(k)

	  do i = 1, m
	    s = 0
	    do j = 1, n
	      s = s + A(i,j)*dconjg(V(k,j))
	    enddo
	    U(k,i) = s*x
	  enddo
	enddo

10	continue

* The vectors corresponding to the zero singular values can be chosen
* arbitrarily.  We choose u_i = delta_ij and use Gram-Schmidt to
* orthonormalize these.

	j = 0
#ifdef DESC
	do i = k, n
#else
	do i = k, 1, -1
#endif
11	  j = j + 1
	  x = 0
	  do r = 1, m
	    if( r .eq. j ) then
	      s = 1
	    else
	      s = 0
	    endif
#ifdef DESC
	    do l = 1, i - 1
#else
	    do l = i + 1, k
#endif
	      s = s - dconjg(U(l,r))*U(l,j)
	    enddo
	    U(i,r) = s
	    x = x + dble(s)**2 + dimag(s)**2
	  enddo
	  if( x .le. eps ) goto 11

	  x = 1/sqrt(x)
	  do r = 1, m
	    U(i,r) = U(i,r)*x
	  enddo
	enddo
	end

