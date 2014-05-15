* PseudoEigensystem.F
* diagonalization of a not quite Hermitian n-by-n matrix
* such that the eigenvalues are nonnegative
* this file is part of FormCalc
* last modified 18 Jun 04 th

* Works only if U and V, the left and right SVD transformation
* matrices, differ by just a phase.  The final transformation matrix
* is U' = diag(sqrt(U V^+)) V.

#include "config.h"


	subroutine PseudoEigensystem(n, A, lda, d, U, ldu)
	implicit none
	integer n, lda, ldu
	double complex A(lda,*), U(ldu,*)
	double precision d(*)

	integer i, j
	double complex V(MAXMAT,MAXMAT), phase

	call SingularValues(n, n, A, lda, d, U, ldu, V, MAXMAT)

	do i = 1, n
	  phase = 0
	  do j = 1, n
	    phase = phase + U(i,j)*dconjg(V(i,j))
	  enddo
	  phase = sqrt(phase)
	  do j = 1, n
	    U(i,j) = phase*V(i,j)
	  enddo
	enddo
	end

