* Li2.f
* the dilogarithm function
* this file is part of FormCalc
* last modified 15 Jun 04 th

* Li2 calls the dilogarithm function of FF, ffzxdl.

	double precision function Li2(x)
	implicit none
	double precision x

	double precision pi
	parameter (pi = 3.1415926535897932384626433832795029D0)

	double complex cli2, dummy
	integer ier, ipi12

	ier = 0
	call ffzxdl(cli2, ipi12, dummy, x, 1, ier)
	Li2 = dble(cli2) + pi**2/12D0*ipi12
	end

