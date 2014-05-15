* ODESolve.F
* ODE solver a la Numerical Recipes
* this file is part of FeynHiggs
* last modified 25 Jul 13 th

#include "externals.h"
#include "types.h"

c#define VERBOSE
c#define WARNINGS


	subroutine ODESolve(nvar, f, ystart, x1, x2, eps, h1)
	implicit none
	integer nvar
	external f
	RealType ystart(*), x1, x2, eps, h1

	RealType safety, pgrow, pshrink, errcon
	parameter (safety = .9D0)
	parameter (pgrow = -.2D0)
	parameter (pshrink = -.25D0)
	parameter (errcon = 1.89D-4)

	RealType A2, B21
	RealType A3, B31, B32
	RealType A4, B41, B42, B43
	RealType A5, B51, B52, B53, B54
	RealType A6, B61, B62, B63, B64, B65
	RealType C1, C3, C4, C6, DC1, DC3, DC4, DC5, DC6
	parameter (A2 = .2D0, B21 = .2D0)
	parameter (A3 = .3D0, B31 = 3/40D0, B32 = 9/40D0)
	parameter (A4 = .6D0, B41 = .3D0, B42 = -.9D0, B43 = 1.2D0)
	parameter (A5 = 1, B51 = -11/54D0, B52 = 2.5D0, B53 = -70/27D0,
     &    B54 = 35/27D0)
	parameter (A6 = .875D0, B61 = 1631/55296D0, B62 = 175/512D0,
     &    B63 = 575/13824D0, B64 = 44275/110592D0, B65 = 253/4096D0)
	parameter (C1 = 37/378D0, DC1 = C1 - 2825/27648D0)
	parameter (C3 = 250/621D0, DC3 = C3 - 18575/48384D0)
	parameter (C4 = 125/594D0, DC4 = C4 - 13525/55296D0)
	parameter (DC5 = -277/14336D0)
	parameter (C6 = 512/1771D0, DC6 = C6 - .25D0)

	integer nvarmax, maxstep
	parameter (nvarmax = 4, maxstep = 10000)

	RealType x0, y0(nvarmax), y(nvarmax), dydx(nvarmax,6)
	RealType h, h0, yerr, yscal, errmax
	integer i, step

	h = sign(h1, x2 - x1)
	x0 = x1
	do i = 1, nvar
	  y0(i) = ystart(i)
	enddo

	do step = 1, maxstep
#ifdef VERBOSE
	  print *, "step,stepsize,x ", step, h, x0
	  print *, "y =", (y0(i), i = 1,nvar)
#endif
	  call f(x0, y0, dydx(1,1))
	  h0 = h
	  if( (x0 + h - x2)*(x0 + h - x1) .gt. 0 ) h = x2 - x0

* determine Cash-Karp step
1	  do i = 1, nvar
	    y(i) = y0(i) + B21*dydx(i,1)*h
	  enddo
	  call f(x0 + A2*h, y, dydx(1,2))

	  do i = 1, nvar
	    y(i) = y0(i) + (B31*dydx(i,1) + B32*dydx(i,2))*h
	  enddo
	  call f(x0 + A3*h, y, dydx(1,3))

	  do i = 1, nvar
	    y(i) = y0(i) + (B41*dydx(i,1) + B42*dydx(i,2) +
     &                      B43*dydx(i,3))*h
	  enddo
	  call f(x0 + A4*h, y, dydx(1,4))

	  do i = 1, nvar
	    y(i) = y0(i) + (B51*dydx(i,1) + B52*dydx(i,2) +
     &                      B53*dydx(i,3) + B54*dydx(i,4))*h
	  enddo
	  call f(x0 + A5*h, y, dydx(1,5))

	  do i = 1, nvar
	    y(i) = y0(i) + (B61*dydx(i,1) + B62*dydx(i,2) +
     &                      B63*dydx(i,3) + B64*dydx(i,4) +
     &                      B65*dydx(i,5))*h
	  enddo
	  call f(x0 + A6*h, y, dydx(1,6))

	  errmax = 0
	  do i = 1, nvar
	    yscal = abs(y0(i)) + abs(h0*dydx(i,1)) + 1D-30
	    yerr = h*(DC1*dydx(i,1) + DC3*dydx(i,3) +
     &                DC4*dydx(i,4) + DC5*dydx(i,5) + DC6*dydx(i,6))
	    errmax = max(errmax, abs(yerr/yscal))
	  enddo

	  errmax = errmax/eps
	  if( errmax .gt. 1 ) then
	    h = sign(max(abs(safety*h*errmax**pshrink), .1D0*abs(h)), h)
	    if( x0 + h .eq. x0 ) then
#ifdef WARNINGS
	      print *, "stepsize underflow in ODESolve"
#endif
	    else
	      goto 1
	    endif
	  endif

	  x0 = x0 + h
	  do i = 1, nvar
	    y0(i) = y0(i) + h*(C1*dydx(i,1) + C3*dydx(i,3) +
     &                         C4*dydx(i,4) + C6*dydx(i,6))
	  enddo

	  if( errmax .gt. errcon ) then
	    h = safety*h*errmax**pgrow
	  else
	    h = 5*h
	  endif

	  if( (x0 - x2)*(x2 - x1) .ge. 0 ) goto 9
	enddo

#ifdef WARNINGS
	print *, "RGE solution did not converge"
#endif

9	do i = 1, nvar
	  ystart(i) = y0(i)
	enddo
	end

