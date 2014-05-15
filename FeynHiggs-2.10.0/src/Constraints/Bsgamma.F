* Bsgamma.F
* computes the rate for B -> X_s gamma
* this file is part of FeynHiggs
* last modified 3 Dec 13 th

#include "externals.h"
#include "types.h"
#include "debug.h"


	subroutine Bsgamma(bsgMSSM, bsgSM)
	implicit none
	RealType bsgMSSM, bsgSM

#include "FH.h"
#include "looptools.h"

	integer Gen4, Cha4, Neu4, All4, Ind1, Ind2, nfdummy
	ComplexType dup1, dup2
	ComplexType tmp1, tmp2, tmp3, tmp4, tmp5
	RealType eta, asMW, asMB

	ComplexType C7L, C7R, C8L, C8R, C7LSM, C8LSM
	ComplexType C7LHp, C7LCha, C7LNeu, C7LGlu
	ComplexType                C7RNeu, C7RGlu
	ComplexType C8LHp, C8LCha, C8LNeu, C8LGlu
	ComplexType                C8RNeu, C8RGlu
	RealType C2, h2, h7, h8, pre

	RealType k77, k771, k78, k88, k22, k27, k28
c	parameter (k77 = 75.67D-2 + 3.47D-2)
	parameter (k77 = 75.67D-2)
	parameter (k771 = 3.47D-2)
	parameter (k22 = .23D-2)
	parameter (k88 = 8.47D-2)
	parameter (k27 = -14.77D-2)
	parameter (k78 = 9.45D-2)
	parameter (k28 = -.04D-2)

	RealType phasespace, z
	phasespace(z) = 1 + z*(z*((8 - z)*z - 12*log(z)) - 8)

	if( abs(CKM(2,3)) .lt. 1D-15 .or. MC .lt. 1D-15 ) then
	  bsgMSSM = 0
	  bsgSM = 0
	  return
	endif

c	call AlphaS(asMT, nfdummy, MT)
c	MTr2 = MT2*(asMW/asMT)**(24/23D0)*(1 - 8/(3*pi)*asMT)
c	MTr = sqrt(MTr2)

	call AlphaS(asMW, nfdummy, MW)
	call AlphaS(asMB, nfdummy, MB)
	eta = asMW/asMB

#include "bsgamma.h"
#include "bsglu.h"

	C7L = C7LSM + C7LHp + C7LCha + C7LNeu + C7LGlu
	C7R =                          C7RNeu + C7RGlu

	C8L = C8LSM + C8LHp + C8LCha + C8LNeu + C8LGlu
	C8R =                          C8RNeu + C8RGlu

	if( debuglevel .ge. 2 ) then
	  DCONST "bsgamma Wilson coefficients before running" ENDL
	  DCONST "C7LSM =", C7LSM	ENDL
	  DCONST "C7L   =", C7L		ENDL
	  DCONST "C7R   =", C7R		ENDL
	  DCONST "C8LSM =", C8LSM	ENDL
	  DCONST "C8L   =", C8L		ENDL
	  DCONST "C8R   =", C8R		ENDL
	endif

	h2 = 626126/272277D0*eta**(14/23D0) -
     &    56281/51730D0*eta**(16/23D0) -
     &    3/7D0*eta**(6/23D0) -
     &    1/14D0*eta**(-12/23D0) -
     &    .6494D0*eta**(.4086D0) -
     &    .0380D0*eta**(-.4230D0) -
     &    .0186D0*eta**(-.8994D0) -
     &    .0057D0*eta**(.1456D0)
	h7 = eta**(16/23D0)
	h8 = 8/3D0*(eta**(14/23D0) - eta**(16/23D0))

	C7L = h7*C7L + h8*C8L + h2
	C7R = h7*C7R + h8*C8R

	C7LSM = h7*C7LSM + h8*C8LSM + h2

	h2 = 313063/363036D0*eta**(14/23D0) -
     &    .9135D0*eta**(.4086D0) +
     &    .0873D0*eta**(-.4230D0) -
     &    .0571D0*eta**(-.8994D0) +
     &    .0209D0*eta**(.1456D0)
	h8 = eta**(14/23D0)

	C8L = h8*C8L + h2
	C8R = h8*C8R

	C8LSM = h8*C8LSM + h2

	C2 = .5D0*(eta**(-12/23D0) + eta**(6/23D0))

	if( debuglevel .ge. 2 ) then
	  DCONST "bsgamma Wilson coefficients after running" ENDL
	  DCONST "C2    =", C2		ENDL
	  DCONST "C7LSM =", C7LSM	ENDL
	  DCONST "C7L   =", C7L		ENDL
	  DCONST "C7R   =", C7R		ENDL
	  DCONST "C8LSM =", C8LSM	ENDL
	  DCONST "C8L   =", C8L		ENDL
	  DCONST "C8R   =", C8R		ENDL
	endif

	pre = .105D0*6/pi*Alfa0/phasespace(MC2/MB2)*
     &    abs(CKMC(3,2)*CKM(3,3)/CKM(2,3))**2

	bsgSM = pre*(
     &    k77*Sq(C7LSM) +
     &    k88*Sq(C8LSM) +
     &    k78*Re(C7LSM*Conjugate(C8LSM)) +
     &    C2*(k22*C2 + k27*Re(C7LSM) + k28*Re(C8LSM)) )

	bsgMSSM = pre*(
     &    k77*(Sq(C7L) + Sq(C7R)) +
     &    k88*(Sq(C8L) + Sq(C8R)) +
     &    k78*(Re(C7L*Conjugate(C8L)) + Re(C7R*Conjugate(C8R))) +
     &    C2*(k22*C2 + k27*Re(C7L) + k28*Re(C8L)) )
	end

