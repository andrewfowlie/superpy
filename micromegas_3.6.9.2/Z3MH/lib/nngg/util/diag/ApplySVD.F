* ApplySVD.F
* transform the n-by-n matrix A by the left and right transformation
* matrices U and V, i.e. form U^* A V^+.
* The original matrix A is overwritten with the result.
* this file is part of FormCalc
* last modified 18 Jun 04 th

#include "config.h"


	subroutine ApplySVD(n, A, lda, U, ldu, V, ldv)
	implicit none
	integer n, lda, ldu, ldv
	double complex A(lda,*), U(ldu,*), V(ldv,*)

	integer i, j, k
	double complex AVdag(MAXMAT,MAXMAT), s

	do j = 1, n
	  do i = 1, n
	    s = 0
	    do k = 1, n
	      s = s + A(i,k)*dconjg(V(j,k))
	    enddo
	    AVdag(i,j) = s
	  enddo
	enddo

	do j = 1, n
	  do i = 1, n
	    s = 0
	    do k = 1, n
	      s = s + dconjg(U(i,k))*AVdag(k,j)
	    enddo
	    A(i,j) = s
	  enddo
	enddo
	end

