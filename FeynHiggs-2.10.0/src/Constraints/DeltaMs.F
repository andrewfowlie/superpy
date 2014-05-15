* DeltaMs.F
* computes the rate for Delta M_{B_s}
* the calculation follows closely hep-ph/0102316
* this file is part of FeynHiggs
* last modified 10 May 12 th

#include "externals.h"
#include "types.h"
#include "debug.h"


	subroutine DeltaMs(deltaMsMSSM, deltaMsSM)
	implicit none
	RealType deltaMsMSSM, deltaMsSM

#include "FH.h"
#include "looptools.h"

	ComplexType CVLL, CVLLSM, CVLLHp, CVLLCha, CVLLNeu, CVLLGlu
	ComplexType CVRR, CVRRNeu, CVRRGlu
	ComplexType CLR1, CLR1Neu, CLR1Glu
	ComplexType CLR2, CLR2Neu, CLR2Glu
	ComplexType CSLL1, CSLL1Neu, CSLL1Glu
	ComplexType CSLL2, CSLL2Neu, CSLL2Glu
	ComplexType CSRR1, CSRR1Neu, CSRR1Glu
	ComplexType CSRR2, CSRR2Neu, CSRR2Glu

	integer Cha5, Cha6, Neu5, Neu6
	integer All5, All6, Ind1, Ind2, Ind3, Ind4
	ComplexType dup1, dup2, dup3, dup4, dup5
	ComplexType dup6, dup7, dup8, tmp1

	integer nfdummy
	RealType asMB, asMT, eta5, nlo, et1, et2, eff
	RealType etaVLL, etaLR(2,2), etaSLL(2,2)
	RealType BVLL, BLR1, BLR2, BSLL1, BSLL2
	RealType PVLL, PLR1, PLR2, PSLL1, PSLL2

	RealType MBs, fBs, hbarc
	parameter (MBs = 5367.5D-3)
	parameter (fBs = 240D-3)
* hbar c in ps*GeV
	parameter (hbarc = 6.58211899D-22*1D12*1D-3)

	if( abs(CKM(2,3)) .lt. 1D-14 ) then
	  deltaMsMSSM = 0
	  deltaMsSM = 0
	  return
	endif

* the implementation follows closely hep-ph/0102316

	call AlphaS(asMB, nfdummy, MB + 1D-14)
	call AlphaS(asMT, nfdummy, Mf(tT,3) - 1D-14)
	eta5 = asMT/asMB

#if 1
* NLO running
	nlo = asMB/(4*pi)

	etaVLL = (1 + nlo*1.6273D0*(1 - eta5))*eta5**(6/23D0)

	et1 = eta5**(3/23D0)
	et2 = eta5**(-24/23D0)

	etaLR(1,1) =
     &    (1 + nlo*(-2.0994D0 + 1.1744D0*eta5))*et1 +
     &    nlo*.9250D0*et2
	etaLR(2,1) =
     &    (2/3D0 + nlo*(-11.7329D0 + .7829D0*eta5))*et1 +
     &    (-2/3D0 + nlo*(-5.3048D0 + 16.2548D0*eta5))*et2
	etaLR(1,2) = nlo*1.3875D0*(eta5*et1 - et2)
	etaLR(2,2) =
     &    nlo*.9250D0*eta5*et1 +
     &    (1 + nlo*(7.9572D0 - 8.8822D0*eta5))*et2

	et1 = eta5**(-.6315D0)
	et2 = eta5**.7184D0

	etaSLL(1,1) =
     &    (1.0153D0 + nlo*(4.8177D0 - 5.2272D0*eta5))*et1 +
     &    (-.0153 + nlo*(.3371D0 + .0724D0*eta5))*et2
	etaSLL(2,1) =
     &    (-.0081D0 + nlo*(.0531D0 + .0415D0*eta5))*et1 +
     &    (.0081D0 + nlo*(-.0566D0 - .0380D0*eta5))*et2
	etaSLL(1,2) =
     &    (1.9325D0 + nlo*(9.1696D0 - 38.8778D0*eta5))*et1 +
     &    (-1.9325D0 + nlo*(42.5021D0 - 12.7939D0*eta5))*et2
	etaSLL(2,2) =
     &    (-.0153 + nlo*(.1011D0 + .3083D0*eta5))*et1 +
     &    (1.0153D0 + nlo*(-7.1314D0 + 6.7219D0*eta5))*et2
#else
* LO running
	etaVLL = eta5**(6/23D0)

	et1 = eta5**(3/23D0)
	et2 = eta5**(-24/23D0)

	etaLR(1,1) = et1
	etaLR(2,1) = 2/3D0*(et1 - et2)
	etaLR(1,2) = 0
	etaLR(2,2) = et2

	et1 = eta5**(-.6315D0)
	et2 = eta5**.7184D0

	etaSLL(1,1) = 1.0153D0*et1 - .0153*et2
	etaSLL(2,1) = .0081D0*(et2 - et1)
	etaSLL(1,2) = 1.9325D0*(et1 - et2)
	etaSLL(2,2) = 1.0153D0*et2 - .0153*et1
#endif

	eff = 1.44D0*(4.4D0/(MB + MS))**2
	BVLL = .61D0*eff
	BLR1 = .96D0*eff
	BLR2 = 1.30D0*eff
	BSLL1 = .76D0*eff
	BSLL2 = .51D0*eff

	PVLL = etaVLL*BVLL
	PLR1 = -1/2D0*etaLR(1,1)*BLR1 + 3/4D0*etaLR(2,1)*BLR2
	PLR2 = -1/2D0*etaLR(1,2)*BLR1 + 3/4D0*etaLR(2,2)*BLR2
	PSLL1 = -5/8D0*etaSLL(1,1)*BSLL1 - 3/2D0*etaSLL(2,1)*BSLL2
	PSLL2 = -5/8D0*etaSLL(2,1)*BSLL1 - 3/2D0*etaSLL(2,2)*BSLL2

	eff = Mf2(tT,3)/MW2

#include "deltaMs.h"

	CVLL = CVLLSM + CVLLHp + CVLLCha + CVLLNeu + CVLLGlu
	CVRR = CVRRNeu + CVRRGlu
	CLR1 = CLR1Neu + CLR1Glu
	CLR2 = CLR1Neu + CLR1Glu
	CSLL1 = CSLL1Neu + CSLL1Glu
	CSLL2 = CSLL2Neu + CSLL2Glu
	CSRR1 = CSRR1Neu + CSRR1Glu
	CSRR2 = CSRR2Neu + CSRR2Glu

	if( debuglevel .ge. 2 ) then
	  DCONST "deltaMs Wilson coefficients before running" ENDL
	  DCONST "CVLLSM =", CVLLSM	ENDL
	  DCONST "CVLL =", CVLL		ENDL
	  DCONST "CVRR =", CVRR		ENDL
	  DCONST "CLR1 =", CLR1		ENDL
	  DCONST "CLR2 =", CLR2		ENDL
	  DCONST "CSLL1 =", CSLL1	ENDL
	  DCONST "CSLL2 =", CSLL2	ENDL
	  DCONST "CSRR1 =", CSRR1	ENDL
	  DCONST "CSRR2 =", CSRR2	ENDL
	endif

	eff = 1/6D0/hbarc*MBs*fBs**2

	deltaMsSM = eff*abs(PVLL*CVLLSM)

	deltaMsMSSM = eff*abs(
     &    PVLL*(CVLL + CVRR) +
     &    PLR1*CLR1 + PLR2*CLR2 +
     &    PSLL1*(CSLL1 + CSRR1) + PSLL2*(CSLL2 + CSRR2) )
	end

