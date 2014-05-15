* VeS.F
* computes Vector.eps.Spinor
* this file is part of FormCalc
* last modified 15 Jun 04 th


	integer function VeS(a, b)
	implicit none
	integer a, b

	double complex vec(2,2, 1)
	common /vectors/ vec

	double complex spi(2, 2)
	equivalence (spi, vec)

#define A(x,y) vec(x,y, a)
#define B(x) spi(x, b)
#define accu 1

	spi(1, accu) = A(1,1)*B(2) - A(1,2)*B(1)
	spi(2, accu) = A(2,1)*B(2) - A(2,2)*B(1)
	VeS = accu
	end

