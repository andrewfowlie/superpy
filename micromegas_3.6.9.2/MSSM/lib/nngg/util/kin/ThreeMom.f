* ThreeMom.f
* this file is part of FormCalc
* last modified 15 Jun 04 th

* ThreeMom computes the length of \vec p_b in the frame in which
* \vec p_a + \vec p_b vanishes.  With sqrtS = sqrt((p_a + p_b)^2),
*
* ThreeMom(sqrtS, ma, mb) = sqrt(lambda(sqrtS, ma**2, mb**2))/(2*sqrtS),
*
* where lambda is the Kallen function.


	double precision function ThreeMom(sqrtS, ma, mb)
	implicit none
	double precision sqrtS, ma, mb

	ThreeMom = sqrt(.25D0*(sqrtS - (ma**2 - mb**2)/sqrtS)**2 -
     &    mb**2)
	end

