* LU.F
* Solution of the linear system A.x = B by LU decomposition
* with partial pivoting
* this file is part of FeynHiggs
* last modified 30 Nov 11 th

* Author: Michael Rauch, 7 Dec 2004
* Reference: Folkmar Bornemann, lecture notes to
* Numerische Mathematik 1, Technical University, Munich, Germany

#include "externals.h"
#include "types.h"
#include "Util.h"

* note: we take A transposed everywhere here for cache efficiency
#undef Element
#define Element(A,i,j) A(j,i)

************************************************************************
* LUDecompose computes the LU decomposition of the n-by-n matrix A
* by Gaussian Elimination with partial pivoting;
* compact (in situ) storage scheme
* Input:
*   A: n-by-n matrix to LU-decompose
*   n: dimension of A
* Output:
*   A: mangled LU decomposition of A in the form
*     ( y11 y12 ... y1n )
*     ( x21 y22 ... y2n )
*     ( x31 x32 ... y3n )
*     ( ............... )
*     ( xn1 xn2 ... ynn )
*   where 
*     (   1   0 ...   0 )  ( y11 y12 ... y1n )
*     ( x21   1 ...   0 )  (   0 y22 ... y2n )
*     ( x31 x32 ...   0 )  (   0   0 ... y3n )  =  Permutation(A)
*     ( ............... )  ( ............... )
*     ( xn1 xn2 ...   1 )  (   0   0 ... ynn ) 
*   perm: permutation vector

	subroutine LUDecompose(n, A,ldA, perm)
	implicit none
	integer n, ldA, perm(*)
	RealType A(ldA,*)

	integer i, j, k, pj, invperm(MAXDIM)
	RealType tmp, absA, pabsA

	do j = 1, n
	  invperm(j) = j
	enddo

	do j = 1, n
* do U part (minus diagonal one)
	  do i = 2, j - 1
	    tmp = 0
	    do k = 1, i - 1
	      tmp = tmp + Element(A,i,k)*Element(A,k,j)
	    enddo
	    Element(A,i,j) = Element(A,i,j) - tmp
	  enddo

* do L part (plus diagonal from U case)
	  pabsA = -1
	  pj = j
	  do i = j, n
	    tmp = 0
	    do k = 1, j - 1
	      tmp = tmp + Element(A,i,k)*Element(A,k,j)
	    enddo
	    Element(A,i,j) = Element(A,i,j) - tmp

* do partial pivoting ...
* find the pivot
	    absA = abs(Element(A,i,j))
	    if( absA .gt. pabsA ) then
	      pabsA = absA
	      pj = i
	    endif
	  enddo

          perm(invperm(pj)) = j

* exchange rows
	  if( pj .ne. j ) then
	    invperm(pj) = invperm(j)
	    do k = 1, n
	      tmp = Element(A,j,k)
	      Element(A,j,k) = Element(A,pj,k)
	      Element(A,pj,k) = tmp
	    enddo
	  endif

* division by the pivot element
	  if( Element(A,j,j) .ne. 0 ) then
	    tmp = 1/Element(A,j,j)
	    do i = j + 1, n
	      Element(A,i,j) = Element(A,i,j)*tmp
	    enddo
	  endif
	enddo
	end


************************************************************************
* LUSolve computes the x in A.x = b from the LU-decomposed A.
* Input:
*   A: LU-decomposed n-by-n matrix A
*   b: input vector b in A.x = b
*   n: dimension of A
*   p: permutation vector from LU decomposition
* Output:
*   b: solution vector x in A.x = b

	subroutine LUSolve(n, A,ldA, b)
	implicit none
	integer n, ldA
	RealType A(ldA,*), b(*)

	integer i, j
	RealType tmp

* forward substitution L.y = b
	do i = 1, n
	  tmp = 0
	  do j = 1, i - 1
	    tmp = tmp + Element(A,i,j)*b(j)
	  enddo
	  b(i) = b(i) - tmp
	enddo

* backward substitution U.x = y
	do i = n, 1, -1
	  if( Element(A,i,i) .ne. 0 ) then
	    tmp = 0
	    do j = i + 1, n
	      tmp = tmp + Element(A,i,j)*b(j)
	    enddo
	    b(i) = (b(i) - tmp)/Element(A,i,i)
	  endif
	enddo
	end

