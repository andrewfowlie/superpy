* RetrieveOSPara.F
* retrieve the on-shell input parameters from FeynHiggs
* this file is part of FeynHiggs
* last modified 30 Nov 11 th

#include "externals.h"
#include "types.h"
#include "debug.h"


	subroutine FHRetrieveOSPara(error, scalefactor_,
     &    MT_, TB_, MA0_, MHp_,
     &    M3SL_, M3SE_, M3SQ_, M3SU_, M3SD_,
     &    M2SL_, M2SE_, M2SQ_, M2SU_, M2SD_,
     &    M1SL_, M1SE_, M1SQ_, M1SU_, M1SD_,
     &    MUE_,
     &    Atau_, At_, Ab_, Amu_, Ac_, As_, Ae_, Au_, Ad_,
     &    M_1_, M_2_, M_3_)
	implicit none
	integer error
	RealType scalefactor_
	RealType MT_, TB_, MA0_, MHp_
	RealType M3SL_, M3SE_, M3SQ_, M3SU_, M3SD_
	RealType M2SL_, M2SE_, M2SQ_, M2SU_, M2SD_
	RealType M1SL_, M1SE_, M1SQ_, M1SU_, M1SD_
	ComplexType MUE_, M_1_, M_2_, M_3_
	ComplexType Atau_, At_, Ab_, Amu_, Ac_, As_, Ae_, Au_, Ad_

#include "FH.h"
#define __SUBROUTINE__ "FHRetrieveOSPara"

	error = 0
	CheckPara(error)

	scalefactor_ = scalefactor

	MT_ = MT
	TB_ = TB
	if( inputmass .eq. 3 ) then
	  MA0_ = MA0
	  MHp_ = -1
	else
	  MA0_ = -1
	  MHp_ = MHp
	endif

	M1SL_ = MSS(1,1)
	M1SE_ = MSS(2,1)
	M1SQ_ = MSS(3,1)
	M1SU_ = MSS(4,1)
	M1SD_ = MSS(5,1)

	M2SL_ = MSS(1,2)
	M2SE_ = MSS(2,2)
	M2SQ_ = MSS(3,2)
	M2SU_ = MSS(4,2)
	M2SD_ = MSS(5,2)

	M3SL_ = MSS(1,3)
	M3SE_ = MSS(2,3)
	M3SQ_ = MSS(3,3)
	M3SU_ = MSS(4,3)
	M3SD_ = MSS(5,3)

	MUE_ = MUE
	M_1_ = M_1
	M_2_ = M_2
	M_3_ = M_3

	Atau_ = Af(2,3)
	At_ = Af(3,3)
	Ab_ = Af(4,3)

	Amu_ = Af(2,2)
	Ac_ = Af(3,2)
	As_ = Af(4,2)

	Ae_ = Af(2,1)
	Au_ = Af(3,1)
	Ad_ = Af(4,1)
	end

