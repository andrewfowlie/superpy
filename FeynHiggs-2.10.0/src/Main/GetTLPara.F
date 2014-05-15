* GetTLPara.F
* get parameters used in the two-loop neutral Higgs-mass calculation
* this file is part of FeynHiggs
* last modified 10 May 12 th

#include "externals.h"
#include "types.h"
#include "debug.h"


	subroutine FHGetTLPara(error,
     &    MSb_, USb_, MbSL2_, Deltab_)
	implicit none
	integer error
	RealType MSb_(2), MbSL2_
	ComplexType USb_(2,2), Deltab_

#include "FH.h"
#define __SUBROUTINE__ "FHGetTLPara"

	error = 0
	CheckSf(error)
	CheckTL(error)
	if( error .gt. 0 ) return

	MSb_(1) = MSf(1,bM0,3)
	MSb_(2) = MSf(2,bM0,3)

	USb_(1,1) = USf(1,1,bM0,3)
	USb_(2,1) = USf(2,1,bM0,3)
	USb_(1,2) = USf(1,2,bM0,3)
	USb_(2,2) = USf(2,2,bM0,3)

	MbSL2_ = DSS2(2,bM0,3) + DSS2(1,bM0,3) - DSf(1,4)
	Deltab_ = Deltaf(bM0,3)
	end

