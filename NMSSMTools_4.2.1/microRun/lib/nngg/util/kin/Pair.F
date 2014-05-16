* Pair.F
* calculates the scalar product of two four-vectors
* this file is part of FormCalc
* last modified 15 Jun 04 th


* Note: the arguments a, b specify the vector's position in
* the array vec, NOT the vector itself.  The macros k, e, ec, s
* (defined in util.h) convert from 1...LEGS to a, b.

	double complex function Pair(a, b)
	implicit none
	integer a, b

	double complex vec(2,2, 1)
	common /vectors/ vec

#define A(x,y) vec(x,y, a)
#define B(x,y) vec(x,y, b)

	Pair = .5D0*(
     &    dble(A(1,1))*dble(B(2,2)) +
     &    dble(A(2,2))*dble(B(1,1)) -
     &    A(1,2)*B(2,1) - A(2,1)*B(1,2) )
	end

