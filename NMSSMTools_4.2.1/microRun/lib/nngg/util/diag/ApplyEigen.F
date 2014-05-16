* ApplyEigen.F
* transform the n-by-n matrix A by a similarity transformation U,
* i.e. form U A U^+.
* The original matrix A is overwritten with the result.
* this file is part of FormCalc
* last modified 18 Jun 04 th

#include "config.h"


	subroutine ApplyEigen(n, A, lda, U, ldu)
	implicit none
	integer n, lda, ldu
	double complex A(lda,*), U(ldu,*)

	integer i, j, k
	double complex AUdag(MAXMAT,MAXMAT), s

	do j = 1, n
	  do i = 1, n
	    s = 0
	    do k = 1, n
	      s = s + A(i,k)*dconjg(U(j,k))
	    enddo
	    AUdag(i,j) = s
	  enddo
	enddo

	do j = 1, n
	  do i = 1, n
	    s = 0
	    do k = 1, n
	      s = s + U(i,k)*AUdag(k,j)
	    enddo
	    A(i,j) = s
	  enddo
	enddo
	end

