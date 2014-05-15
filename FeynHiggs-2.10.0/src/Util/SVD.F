* SVD.F
* singular value decomposition of an m-by-n matrix
* this file is part of the Diag library
* last modified 30 Nov 11 th

#include "externals.h"
#include "types.h"
#include "Util.h"


************************************************************************
** SVD performs a singular value decomposition.
** Input: m, n, A = m-by-n matrix.
** Output: d = nm-vector of singular values,
** V = nm-by-m left transformation matrix,
** W = nm-by-n right transformation matrix, nm = min(m, n),
** these fulfill diag(d) = V^* A W^+.

	subroutine SVD(m, n, Ao,ldA, d, Vo,ldV, Wo,ldW, sort)
	implicit none
	integer m, n, ldA, ldV, ldW, sort
	ComplexType Ao(ldA,*), Vo(ldV,*), Wo(ldW,*)
	RealType d(*)

	integer nx, nm, p, q, px, qx, j, rev, pi(MAXDIM)
	RealType red, off, thresh
	RealType t, dv, dw, xv, xw, invc
	ComplexType App, Apq, Aqp, Aqq
	ComplexType x, y, sv, sw, tv, tw, f
	ComplexType VW(MAXDIM,MAXDIM,0:2)

* note: for better cache efficiency, the Vx, Wx arrays
* contain the *transpose* of the transformation matrices
	ComplexType V(MAXDIM,MAXDIM)
	ComplexType W(MAXDIM,MAXDIM)
	ComplexType A(MAXDIM,MAXDIM)
	equivalence (VW(1,1,0), V)
	equivalence (VW(1,1,1), W)
	equivalence (VW(1,1,2), A)

	integer sweep
	common /nsweeps/ sweep

	nx = max(m, n)

	if( nx .gt. MAXDIM ) then
	  print *, "Dimension too large"
	  d(1) = -999
	  return
	endif

	do p = 1, nx
	  do q = 1, nx
	    V(q,p) = 0
	    W(q,p) = 0
	    A(q,p) = 0
	  enddo
	  V(p,p) = 1
	  W(p,p) = 1
	enddo

	rev = ibits(m - n, 15, 1)
	if( rev .eq. 1 ) then
	  do p = 1, n
	    do q = 1, m
	      A(p,q) = Ao(q,p)
	    enddo
	  enddo
	else
	  do p = 1, n
	    do q = 1, m
	      A(q,p) = Ao(q,p)
	    enddo
	  enddo
	endif

	red = .01D0/nx**4

	do sweep = 1, 50
	  off = 0
	  do q = 2, nx
	    do p = 1, q - 1
	      off = off + Sq(A(p,q)) + Sq(A(q,p))
	    enddo
	  enddo
	  if( .not. off .gt. EPS ) goto 1

	  thresh = 0
	  if( sweep .lt. 4 ) thresh = off*red

	  do q = 2, nx
	    do p = 1, q - 1
	      px = p
	      qx = q
	      if( Sq(A(p,p)) + Sq(A(q,q)) .lt.
     &            Sq(A(p,q)) + Sq(A(q,p)) ) then
	        px = q
	        qx = p
	      endif

	      App = A(px,p)
	      Aqq = A(qx,q)
	      Apq = A(px,q)
	      Aqp = A(qx,p)
	      off = Sq(Apq) + Sq(Aqp)
	      if( sweep .gt. 4 .and.
     &            off .lt. EPS*(Sq(App) + Sq(Aqq)) ) then
	        A(px,q) = 0
	        A(qx,p) = 0
	      else if( off .gt. thresh ) then
	        xv = Re((App - Aqq)*Conjugate(App + Aqq))
	        xw = Re((Apq - Aqp)*Conjugate(Apq + Aqp))
	        dv = .5D0*(xv + xw)
	        dw = .5D0*(xv - xw)

	        tv = Conjugate(App)*Aqp + Aqq*Conjugate(Apq)
	        tw = Conjugate(App)*Apq + Aqq*Conjugate(Aqp)
c	        t = sqrt(dv**2 + Sq(tv))
	        t = sqrt(dw**2 + Sq(tw))

	        xv = min(abs(dv + t), abs(dw + t))
	        xw = min(abs(dv - t), abs(dw - t))
	        if( xv + xw .gt. DBL_EPS ) then
	          t = sign(t, xv - xw)
	          tv = tv/(dv + t)
	          tw = tw/(dw + t)
	        else
	          tv = 0
	          tw = Apq/App
	        endif

	        invc = sqrt(1 + Sq(tv))
	        sv = tv/invc
	        tv = tv/(invc + 1)

	        invc = sqrt(1 + Sq(tw))
	        sw = tw/invc
	        tw = tw/(invc + 1)

	        do j = 1, nx
	          x = A(j,p)
	          y = A(j,q)
	          A(j,p) = x + Conjugate(sw)*(y - tw*x)
	          A(j,q) = y - sw*(x + Conjugate(tw)*y)
	          x = A(px,j)
	          y = A(qx,j)
	          A(p,j) = x + Conjugate(sv)*(y - tv*x)
	          A(q,j) = y - sv*(x + Conjugate(tv)*y)
	        enddo

	        A(p,p) = invc*(App + Conjugate(sv)*(Aqp - tv*App))
	        A(q,p) = 0
	        A(p,q) = 0
	        A(q,q) = invc*(Aqq - sv*(Apq + Conjugate(tv)*Aqq))

	        do j = 1, nx
	          x = V(j,px)
	          y = V(j,qx)
	          V(j,p) = x + sv*(y - Conjugate(tv)*x)
	          V(j,q) = y - Conjugate(sv)*(x + tv*y)
	        enddo

	        do j = 1, nx
	          x = W(j,p)
	          y = W(j,q)
	          W(j,p) = x + sw*(y - Conjugate(tw)*x)
	          W(j,q) = y - Conjugate(sw)*(x + tw*y)
	        enddo
	        goto 2
	      endif

	      if( p .ne. px ) then
	        do j = 1, nx
	          x = A(p,j)
	          A(p,j) = A(q,j)
	          A(q,j) = x
	        enddo

	        do j = 1, nx
	          x = V(j,p)
	          V(j,p) = V(j,q)
	          V(j,q) = x
	        enddo
	      endif

2	      continue
	    enddo
	  enddo
	enddo

	print *, "Bad convergence in SVD"

1	continue

	nm = min(m, n)

* make the diagonal elements nonnegative

	do p = 1, nm
	  d(p) = abs(A(p,p))
	  if( d(p) .gt. DBL_EPS .and. d(p) .ne. Re(A(p,p)) ) then
	    f = A(p,p)/d(p)
	    do q = 1, nm
	      W(q,p) = W(q,p)*f
	    enddo
	  endif
	enddo

* sort the singular values

	do p = 1, nm
	  pi(p) = p
	enddo

	do p = 1, nm
	  j = p
	  t = d(p)
	  if( sort .ne. 0 ) then
	    do q = p + 1, nm
	      if( sort*(t - d(q)) .gt. 0 ) then
	        j = q
	        t = d(q)
	      endif
	    enddo
	  endif

	  d(j) = d(p)
	  d(p) = t

	  q = pi(j)
	  pi(j) = pi(p)

	  do j = 1, m
	    Vo(p,j) = VW(j,q,rev)
	  enddo
	  do j = 1, n
	    Wo(p,j) = VW(j,q,1-rev)
	  enddo
	enddo
	end

