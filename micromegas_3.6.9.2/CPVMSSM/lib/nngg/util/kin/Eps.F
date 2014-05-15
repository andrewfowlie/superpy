* Eps.F
* calculates -I times the Levi-Civita tensor contracted with
* four four-vectors
* this file is part of FormCalc
* last modified 15 Jun 04 th


* Note: the arguments a, b, c, d specify the vector's position in
* the array vec, NOT the vector itself.  The macros k, e, ec, s
* (defined in util.h) convert from 1...LEGS to a, b, c, d.

	double complex function Eps(a, b, c, d)
	implicit none
	integer a, b, c, d

	double complex vec(2,2, 1)
	common /vectors/ vec

#define A(x,y) vec(x,y, a)
#define B(x,y) vec(x,y, b)
#define C(x,y) vec(x,y, c)
#define D(x,y) vec(x,y, d)

#define diag(i, j) (dble(i(1,1))*dble(j(2,2)) - dble(i(2,2))*dble(j(1,1)))
#define offdiag(i, j) (i(2,1)*j(1,2) - i(1,2)*j(2,1))

	Eps = .25D0*(
     &    diag(A, B)*
     &      offdiag(C, D) +
     &    diag(C, D)*
     &      offdiag(A, B) -
     &    diag(A, C)*
     &      offdiag(B, D) -
     &    diag(B, D)*
     &      offdiag(A, C) +
     &    diag(A, D)*
     &      offdiag(B, C) +
     &    diag(B, C)*
     &      offdiag(A, D) )
	end

