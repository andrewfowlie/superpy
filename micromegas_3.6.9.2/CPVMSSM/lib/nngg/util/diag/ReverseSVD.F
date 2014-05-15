* ReverseSVD.F
* recover the original matrix A from the singular values and
* the left and right transformation matrices U and V by
* forming A = U^T diag(d) V.
* this file is part of FormCalc
* last modified 18 Jun 04 th


	subroutine ReverseSVD(n, A, lda, d, U, ldu, V, ldv)
	implicit none
	integer n, lda, ldu, ldv
	double complex A(lda,*), U(ldu,*), V(ldv,*)
	double precision d(*)

	integer i, j, k
	double complex s

	do j = 1, n
	  do i = 1, n
	    s = 0
	    do k = 1, n
	      s = s + d(k)*U(k,i)*V(k,j)
	    enddo
	    A(i,j) = s
	  enddo
	enddo
	end

