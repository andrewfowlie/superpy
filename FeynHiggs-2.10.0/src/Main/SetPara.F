* SetPara.F
* set the input parameters for FeynHiggs
* this file is part of FeynHiggs
* last modified 30 Nov 11 th

#include "externals.h"
#include "types.h"


	subroutine FHSetPara(error, scalefactor_,
     &    MT_, TB_, MA0_, MHp_,
     &    M3SL_, M3SE_, M3SQ_, M3SU_, M3SD_,
     &    M2SL_, M2SE_, M2SQ_, M2SU_, M2SD_,
     &    M1SL_, M1SE_, M1SQ_, M1SU_, M1SD_,
     &    MUE_,
     &    Atau_, At_, Ab_, Amu_, Ac_, As_, Ae_, Au_, Ad_,
     &    M_1_, M_2_, M_3_,
     &    Qtau_, Qt_, Qb_)
	implicit none
	integer error
	RealType scalefactor_
	RealType MT_, TB_, MA0_, MHp_
	RealType M3SL_, M3SE_, M3SQ_, M3SU_, M3SD_
	RealType M2SL_, M2SE_, M2SQ_, M2SU_, M2SD_
	RealType M1SL_, M1SE_, M1SQ_, M1SU_, M1SD_
	ComplexType MUE_, M_1_, M_2_, M_3_
	ComplexType Atau_, At_, Ab_, Amu_, Ac_, As_, Ae_, Au_, Ad_
	RealType Qtau_, Qt_, Qb_

#include "FH.h"

	if( sm_valid .ne. valid ) call FHSetSMPara(error,
     &    -1D0, -1D0, -1D0,
     &    -1D0, -1D0, -1D0, -1D0, -1D0, -1D0, -1D0, -1D0,
     &    -1D0, -1D0,
     &    -1D0, -1D0, -1D0, -1D0)

	call BeginPara(error, scalefactor_, MT_, TB_, M_3_)

	call HiggsMasses(error, MA0_, MHp_, 0D0)

	call CharginoMasses(M_2_, MUE_)

	call NeutralinoMasses(M_1_)

	call SfermionMasses(2,
     &    M1SL_, M2SL_, M3SL_,
     &    M1SE_, M2SE_, M3SE_,
     &    Ae_, Amu_, Atau_, Qtau_)
	call SfermionMasses(3,
     &    M1SQ_, M2SQ_, M3SQ_,
     &    M1SU_, M2SU_, M3SU_,
     &    Au_, Ac_, At_, Qt_)
	call SfermionMasses(4,
     &    M1SQ_, M2SQ_, M3SQ_,
     &    M1SD_, M2SD_, M3SD_,
     &    Ad_, As_, Ab_, Qb_)

	call EndPara(error)
	end

