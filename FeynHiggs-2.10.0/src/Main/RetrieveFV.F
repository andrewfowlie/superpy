* RetrieveFV.F
* retrieve the flavour-violation parameters
* this file is part of FeynHiggs
* last modified 29 May 12 th

#include "externals.h"
#include "types.h"
#include "debug.h"


	subroutine FHRetrieveLFV(error,
     &    deltaLLL12_, deltaLLL23_, deltaLLL13_,
     &    deltaELR12_, deltaELR23_, deltaELR13_,
     &    deltaERL12_, deltaERL23_, deltaERL13_,
     &    deltaERR12_, deltaERR23_, deltaERR13_)
	implicit none
	integer error
	ComplexType deltaLLL12_, deltaLLL23_, deltaLLL13_
	ComplexType deltaELR12_, deltaELR23_, deltaELR13_
	ComplexType deltaERL12_, deltaERL23_, deltaERL13_
	ComplexType deltaERR12_, deltaERR23_, deltaERR13_

#include "FH.h"
#define __SUBROUTINE__ "FHRetrieveNMFV"

	error = 0
	CheckPara(error)

	deltaLLL12_ = deltaSf_LL(1,2,2)
	deltaLLL23_ = deltaSf_LL(2,3,2)
	deltaLLL13_ = deltaSf_LL(1,3,2)

	deltaELR12_ = deltaSf_LR(1,2,2)
	deltaELR23_ = deltaSf_LR(2,3,2)
	deltaELR13_ = deltaSf_LR(1,3,2)

	deltaERL12_ = deltaSf_RL(1,2,2)
	deltaERL23_ = deltaSf_RL(2,3,2)
	deltaERL13_ = deltaSf_RL(1,3,2)

	deltaERR12_ = deltaSf_RR(1,2,2)
	deltaERR23_ = deltaSf_RR(2,3,2)
	deltaERR13_ = deltaSf_RR(1,3,2)
	end

#undef __SUBROUTINE__


************************************************************************

	subroutine FHRetrieveNMFV(error,
     &    deltaQLL12_, deltaQLL23_, deltaQLL13_,
     &    deltaULR12_, deltaULR23_, deltaULR13_,
     &    deltaURL12_, deltaURL23_, deltaURL13_,
     &    deltaURR12_, deltaURR23_, deltaURR13_,
     &    deltaDLR12_, deltaDLR23_, deltaDLR13_,
     &    deltaDRL12_, deltaDRL23_, deltaDRL13_,
     &    deltaDRR12_, deltaDRR23_, deltaDRR13_)
	implicit none
	integer error
	ComplexType deltaQLL12_, deltaQLL23_, deltaQLL13_
	ComplexType deltaULR12_, deltaULR23_, deltaULR13_
	ComplexType deltaURL12_, deltaURL23_, deltaURL13_
	ComplexType deltaURR12_, deltaURR23_, deltaURR13_
	ComplexType deltaDLR12_, deltaDLR23_, deltaDLR13_
	ComplexType deltaDRL12_, deltaDRL23_, deltaDRL13_
	ComplexType deltaDRR12_, deltaDRR23_, deltaDRR13_

#include "FH.h"
#define __SUBROUTINE__ "FHRetrieveNMFV"

	error = 0
	CheckPara(error)

c	deltaQLL12_ = deltaSf_LL(1,2,3)
c	deltaQLL23_ = deltaSf_LL(2,3,3)
c	deltaQLL13_ = deltaSf_LL(1,3,3)

	deltaULR12_ = deltaSf_LR(1,2,3)
	deltaULR23_ = deltaSf_LR(2,3,3)
	deltaULR13_ = deltaSf_LR(1,3,3)

	deltaURL12_ = deltaSf_RL(1,2,3)
	deltaURL23_ = deltaSf_RL(2,3,3)
	deltaURL13_ = deltaSf_RL(1,3,3)

	deltaURR12_ = deltaSf_RR(1,2,3)
	deltaURR23_ = deltaSf_RR(2,3,3)
	deltaURR13_ = deltaSf_RR(1,3,3)

	deltaQLL12_ = deltaSf_LL(1,2,4)
	deltaQLL23_ = deltaSf_LL(2,3,4)
	deltaQLL13_ = deltaSf_LL(1,3,4)

	deltaDLR12_ = deltaSf_LR(1,2,4)
	deltaDLR23_ = deltaSf_LR(2,3,4)
	deltaDLR13_ = deltaSf_LR(1,3,4)

	deltaDRL12_ = deltaSf_RL(1,2,4)
	deltaDRL23_ = deltaSf_RL(2,3,4)
	deltaDRL13_ = deltaSf_RL(1,3,4)

	deltaDRR12_ = deltaSf_RR(1,2,4)
	deltaDRR23_ = deltaSf_RR(2,3,4)
	deltaDRR13_ = deltaSf_RR(1,3,4)
	end

