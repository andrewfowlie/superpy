* VecSet.F
* construct the CMS representation of a momentum and the
* three corresponding polarization vectors
* this file is part of FormCalc
* last modified 15 Jun 04 th


* i: the index of the momentum
* m, p: mass and three-momentum of the particle
* ex, ey, ez: the unit three-vector of the momentum

	subroutine VecSet(i, m, p, ex, ey, ez)
	implicit none
	integer i
	double precision m, p, ex, ey, ez

	double complex vec(2,2, 8, 0:1)
	common /vectors/ vec

#define K(x,y)     vec(x,y, 1, i)
#define EpsTR(x,y) vec(x,y, 2, i)
#define EpsL(x,y)  vec(x,y, 3, i)
#define EpsTL(x,y) vec(x,y, 4, i)

	double complex spi(2, 16, 0:1)
	equivalence (spi, vec)

#define Spinor(hel, om, x) spi(x, hel*2 + om + 5, i)
#define DottedSpinor(hel, om, x) spi(x, hel*2 + om + 7, i)

* i.e., spi(x,  9, i) = left-handed,  undotted, upper
*       spi(x, 10, i) = left-handed,  undotted, lower
*       spi(x, 11, i) = left-handed,  dotted,   upper
*       spi(x, 12, i) = left-handed,  dotted,   lower
*       spi(x, 13, i) = right-handed, undotted, upper
*       spi(x, 14, i) = right-handed, undotted, lower
*       spi(x, 15, i) = right-handed, dotted,   upper
*       spi(x, 16, i) = right-handed, dotted,   lower

	double precision momspec(4, 1)
	common /momenta/ momspec

	double precision m2, p2, p0, deltap, sump, sinth
	double complex expIfi

	double precision sqrt2
	parameter (sqrt2 = 1.4142135623730950488016887242096981D0)

	m2 = m**2
	if( m2 .lt. 1D-14 ) then
	  p0 = p
	  deltap = 0
	else
	  p2 = p**2
	  p0 = sqrt(p2 + m2)
	  deltap = p0 - p
	  if( deltap .lt. 1D-8*p ) deltap = m2/p*(.5D0 - .125D0*m2/p2)
	endif

	momspec(1, i) = m2
	momspec(2, i) = p0
	momspec(3, i) = p
	momspec(4, i) = deltap

	K(1,1) = p0*(1 + ez) - deltap*ez
	K(2,2) = p0*(1 - ez) + deltap*ez
	K(2,1) = p*dcmplx(ex, ey)
	K(1,2) = dconjg(K(2,1))

	if( m .ne. 0 ) then
	  EpsL(1,1) = (p*(1 + ez) + deltap*ez)/m
	  EpsL(2,2) = (p*(1 - ez) - deltap*ez)/m
	  EpsL(2,1) = p0/m*dcmplx(ex, ey)
	  EpsL(1,2) = dconjg(EpsL(2,1))
	endif

	sinth = sqrt(ex**2 + ey**2)
	if( sinth .lt. 1D-14 ) then
* phi is irrelevant when theta = 0
	  expIfi = 1
	else
	  expIfi = dcmplx(ex, ey)/sinth
	endif

	EpsTR(2,2) = 1/sqrt2*sinth
	EpsTL(2,2) = EpsTR(2,2)
	EpsTR(1,1) = -EpsTR(2,2)
	EpsTL(1,1) = EpsTR(1,1)
	EpsTR(2,1) = 1/sqrt2*(ez - 1)*expIfi
	EpsTL(1,2) = dconjg(EpsTR(2,1))
	EpsTL(2,1) = 1/sqrt2*(ez + 1)*expIfi
	EpsTR(1,2) = dconjg(EpsTL(2,1))

* this is E^(I fi) cos(th/2) = 1/sqrt2*sqrt(1 + ez)*expIfi:
	expIfi = 1/sqrt2*sqrt(1 + ez)*expIfi
* this is sin(th/2):
	sinth = 1/sqrt2*sqrt(1 - ez)
	sump = sqrt(p0 + p)
	deltap = sqrt(deltap)

	Spinor(-1, 6, 1) = deltap*sinth
	DottedSpinor(-1, 6, 1) = Spinor(-1, 6, 1)
	Spinor(-1, 6, 2) = -deltap*expIfi
	DottedSpinor(-1, 6, 2) = dconjg(Spinor(-1, 6, 2))

	Spinor(-1, 7, 1) = sump*sinth
	DottedSpinor(-1, 7, 1) = Spinor(-1, 7, 1)
	Spinor(-1, 7, 2) = -sump*expIfi
	DottedSpinor(-1, 7, 2) = dconjg(Spinor(-1, 7, 2))

	DottedSpinor(1, 6, 1) = sump*expIfi
	Spinor(1, 6, 1) = dconjg(DottedSpinor(1, 6, 1))
	DottedSpinor(1, 6, 2) = sump*sinth
	Spinor(1, 6, 2) = DottedSpinor(1, 6, 2)

	DottedSpinor(1, 7, 1) = deltap*expIfi
	Spinor(1, 7, 1) = dconjg(DottedSpinor(1, 7, 1))
	DottedSpinor(1, 7, 2) = deltap*sinth
	Spinor(1, 7, 2) = DottedSpinor(1, 7, 2)
	end

