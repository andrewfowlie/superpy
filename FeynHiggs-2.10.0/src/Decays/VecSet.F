* VecSet.F
* construct the CMS representation of a momentum and the
* three corresponding polarization vectors
* this file is part of FeynHiggs
* last modified 21 Sep 12 th

#include "externals.h"
#include "types.h"


* i: the index of the momentum
* m, p: mass and three-momentum of the particle
* ex,ey,ez: the unit three-vector of the momentum

	subroutine VecSet(i, m, p, ex,ey,ez)
	implicit none
	integer i
	RealType m, p, ex,ey,ez

#include "Decay.h"

	ComplexType vec_(2,2, nvec, LEGS)
	equivalence (vec, vec_)

#define K(x,y)     vec_(x,y, 1, i)
#define EpsTR(x,y) vec_(x,y, 2, i)
#define EpsL(x,y)  vec_(x,y, 3, i)
#define EpsTL(x,y) vec_(x,y, 4, i)

	ComplexType spi_(2, 6:7, nvec, 1)
	equivalence (vec, spi_)

#define Spi(hel,om,x) spi_(x, om, hel+6, i)
#define DotSpi(hel,om,x) spi_(x, om, hel+7, i)

* i.e., spi(x, 6, 5, i) = left-handed,  undotted, upper
*       spi(x, 7, 5, i) = left-handed,  undotted, lower
*       spi(x, 6, 6, i) = left-handed,  dotted,   upper
*       spi(x, 7, 6, i) = left-handed,  dotted,   lower
*       spi(x, 6, 7, i) = right-handed, undotted, upper
*       spi(x, 7, 7, i) = right-handed, undotted, lower
*       spi(x, 6, 8, i) = right-handed, dotted,   upper
*       spi(x, 7, 8, i) = right-handed, dotted,   lower

	RealType p2, p0, deltap, sump
	RealType sinth, onePez, oneMez
	ComplexType expIphi

	mass(1,i) = m
	mass(2,i) = m**2
	if( mass(2,i) .lt. 1D-14 ) then
	  p0 = p
	  deltap = 0
	else
	  p2 = p**2
	  p0 = sqrt(p2 + mass(2,i))
	  deltap = mass(2,i)/(p0 + p)
	endif

	sinth = ex**2 + ey**2
	onePez = 1 + ez
	if( onePez .lt. .5D0 ) onePez = sinth/(1 - ez)
	oneMez = 1 - ez
	if( oneMez .lt. .5D0 ) oneMez = sinth/(1 + ez)

	if( sinth .lt. 1D-14 ) then
* phi is irrelevant when theta = 0
	  expIphi = 1
	else
	  sinth = sqrt(sinth)
	  expIphi = ToComplex(ex, ey)/sinth
	endif

	K(1,1) = p0*onePez - deltap*ez
	K(2,2) = p0*oneMez + deltap*ez
	K(2,1) = p*ToComplex(ex, ey)
	K(1,2) = Conjugate(K(2,1))

	if( m .ne. 0 ) then
	  EpsL(1,1) = (p*onePez + deltap*ez)/m
	  EpsL(2,2) = (p*oneMez - deltap*ez)/m
	  EpsL(2,1) = p0/m*ToComplex(ex, ey)
	  EpsL(1,2) = Conjugate(EpsL(2,1))
	endif

	EpsTR(2,2) = 1/sqrt2*sinth
	EpsTL(2,2) = EpsTR(2,2)
	EpsTR(1,1) = -EpsTR(2,2)
	EpsTL(1,1) = EpsTR(1,1)
	EpsTR(2,1) = -1/sqrt2*oneMez*expIphi
	EpsTL(1,2) = Conjugate(EpsTR(2,1))
	EpsTL(2,1) = 1/sqrt2*onePez*expIphi
	EpsTR(1,2) = Conjugate(EpsTL(2,1))

* this is E^(I phi) cos(th/2) = 1/sqrt2*sqrt(1 + ez)*expIphi:
	expIphi = 1/sqrt2*sqrt(onePez)*expIphi
* this is sin(th/2):
	sinth = 1/sqrt2*sqrt(oneMez)
	sump = sqrt(p0 + p)
	deltap = sqrt(deltap)

	Spi(-1, 6, 1) = deltap*sinth
	DotSpi(-1, 6, 1) = Spi(-1, 6, 1)
	Spi(-1, 6, 2) = -deltap*expIphi
	DotSpi(-1, 6, 2) = Conjugate(Spi(-1, 6, 2))

	Spi(-1, 7, 1) = sump*sinth
	DotSpi(-1, 7, 1) = Spi(-1, 7, 1)
	Spi(-1, 7, 2) = -sump*expIphi
	DotSpi(-1, 7, 2) = Conjugate(Spi(-1, 7, 2))

	DotSpi(1, 6, 1) = sump*expIphi
	Spi(1, 6, 1) = Conjugate(DotSpi(1, 6, 1))
	DotSpi(1, 6, 2) = sump*sinth
	Spi(1, 6, 2) = DotSpi(1, 6, 2)

	DotSpi(1, 7, 1) = deltap*expIphi
	Spi(1, 7, 1) = Conjugate(DotSpi(1, 7, 1))
	DotSpi(1, 7, 2) = deltap*sinth
	Spi(1, 7, 2) = DotSpi(1, 7, 2)
	end

