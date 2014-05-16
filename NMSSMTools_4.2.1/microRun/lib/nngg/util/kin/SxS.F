* SxS.F
* computes Spinor.Spinor
* this file is part of FormCalc
* last modified 15 Jun 04 th


	double complex function SxS(a, b)
	implicit none
	integer a, b

	double complex spi(2, 1)
	common /vectors/ spi

#define A(x) spi(x, a)
#define B(x) spi(x, b)

	SxS = A(1)*B(1) + A(2)*B(2)
	end

