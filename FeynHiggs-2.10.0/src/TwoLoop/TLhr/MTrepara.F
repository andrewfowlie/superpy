* MTrepara.F
* terms correcting the use of MT-DRbar at 1L
* arXiv:1103.1335 Eqs. (95)-(102) and App. C
* this file is part of FeynHiggs
* last modified 30 Nov 11 th

#include "externals.h"
#include "types.h"


	subroutine MTrepara(imass, dMT, dPhi)
	implicit none
	integer imass
	RealType dMT, dPhi(6)

#include "TLhr.h"
#define __SUBROUTINE__ "MTrepara"

	RealType MSQ2, Xt2, ReAtCXt, dl22, dcH, dcA, dlA, G, c
	RealType dsig11, dsig12, dsig22, dsig1A, dsig2A, dsigAA

	MSQ2 = Re(MSS2(3,3,3))
	Xt2 = Sq(Xt)
	ReAtCXt = Re(AtC*Xt)

#include "MTrepara.h"

	if( imass .eq. 3 ) then
	  G = -Im(MUEXt)**2*dcA
	else
	  G = .5D0*Sq(MUE)*dcH
	endif

	c = 3/(4*pi)*Alfa2L/(MW2*SW2*SB2) * 2*MTy
	dsig11 = c*(G + Re(MUEXt)**2*dcA)
	dsig12 = c*(-G/TB - Re(MUEXt)*(ReAtCXt*dcA - dlA))
	dsig22 = c*(-dl22 + G/TB2 + ReAtCXt*(ReAtCXt*dcA - 2*dlA))
	c = c/SB
	dsig1A = c*(-Im(MUEXt**2)*dcA)
	dsig2A = c*(Im(MUEXt)*(ReAtCXt*dcA - dlA))
	c = c/SB
	dsigAA = c*(G + Im(MUEXt)**2*dcA)

	dPhi(se11) = dPhi(se11) + dMT*dsig11
	dPhi(se12) = dPhi(se12) + dMT*dsig12
	dPhi(se22) = dPhi(se22) + dMT*dsig22
	dPhi(se1A) = dPhi(se1A) + dMT*dsig1A
	dPhi(se2A) = dPhi(se2A) + dMT*dsig2A
	dPhi(seAA) = dPhi(seAA) + dMT*dsigAA

	if( debuglevel .ge. 4 ) then
	  DSELF "dMT shifts:" ENDL
	  DSELF "dsig11 =", dsig11 ENDL
	  DSELF "dsig12 =", dsig12 ENDL
	  DSELF "dsig22 =", dsig22 ENDL
	  DSELF "dsig1A =", dsig1A ENDL
	  DSELF "dsig2A =", dsig2A ENDL
	  DSELF "dsigAA =", dsigAA ENDL
	endif
	end

