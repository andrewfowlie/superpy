* SetFV.F
* set the flavour-violation parameters
* this file is part of FeynHiggs
* last modified 25 Apr 13 th

#include "externals.h"
#include "types.h"
#include "debug.h"


	subroutine FHSetLFV(error,
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
#include "FHRecord.h"
#define __SUBROUTINE__ "FHSetLFV"

	error = 0
	CheckPara(error)
	if( error .gt. 0 ) return

	if( debuglevel .ge. 2 ) then
	  DPARA "============================"	ENDL
	  DPARA "FHSetLFV input:"		ENDL
	  DPARA "----------------------------"	ENDL

	  DTAGrez(ideltaLLL12, deltaLLL12_)
	  DTAGimz(ideltaLLL12, deltaLLL12_)
	  DTAGrez(ideltaLLL23, deltaLLL23_)
	  DTAGimz(ideltaLLL23, deltaLLL23_)
	  DTAGrez(ideltaLLL13, deltaLLL13_)
	  DTAGimz(ideltaLLL13, deltaLLL13_)

	  DTAGrez(ideltaELR12, deltaELR12_)
	  DTAGimz(ideltaELR12, deltaELR12_)
	  DTAGrez(ideltaELR23, deltaELR23_)
	  DTAGimz(ideltaELR23, deltaELR23_)
	  DTAGrez(ideltaELR13, deltaELR13_)
	  DTAGimz(ideltaELR13, deltaELR13_)

	  DTAGrez(ideltaERL12, deltaERL12_)
	  DTAGimz(ideltaERL12, deltaERL12_)
	  DTAGrez(ideltaERL23, deltaERL23_)
	  DTAGimz(ideltaERL23, deltaERL23_)
	  DTAGrez(ideltaERL13, deltaERL13_)
	  DTAGimz(ideltaERL13, deltaERL13_)

	  DTAGrez(ideltaERR12, deltaERR12_)
	  DTAGimz(ideltaERR12, deltaERR12_)
	  DTAGrez(ideltaERR23, deltaERR23_)
	  DTAGimz(ideltaERR23, deltaERR23_)
	  DTAGrez(ideltaERR13, deltaERR13_)
	  DTAGimz(ideltaERR13, deltaERR13_)
	endif

	sf_valid = 0
	tl_valid = 0
	higgs_valid = 0
	coup_valid = 0

	deltaSf_LL(1,2,2) = deltaLLL12_
	deltaSf_LL(2,3,2) = deltaLLL23_
	deltaSf_LL(1,3,2) = deltaLLL13_

	deltaSf_LR(1,2,2) = deltaELR12_
	deltaSf_LR(2,3,2) = deltaELR23_
	deltaSf_LR(1,3,2) = deltaELR13_

	deltaSf_RL(1,2,2) = deltaERL12_
	deltaSf_RL(2,3,2) = deltaERL23_
	deltaSf_RL(1,3,2) = deltaERL13_

	deltaSf_RR(1,2,2) = deltaERR12_
	deltaSf_RR(2,3,2) = deltaERR23_
	deltaSf_RR(1,3,2) = deltaERR13_

	call IniFV(2)
	end

#undef __SUBROUTINE__


************************************************************************

	subroutine FHSetNMFV(error,
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
#include "FHRecord.h"
#define __SUBROUTINE__ "FHSetNMFV"

	error = 0
	CheckPara(error)
	if( error .gt. 0 ) return

	if( debuglevel .ge. 2 ) then
	  DPARA "============================"	ENDL
	  DPARA "FHSetNMFV input:"		ENDL
	  DPARA "----------------------------"	ENDL

	  DTAGrez(ideltaQLL12, deltaQLL12_)
	  DTAGimz(ideltaQLL12, deltaQLL12_)
	  DTAGrez(ideltaQLL23, deltaQLL23_)
	  DTAGimz(ideltaQLL23, deltaQLL23_)
	  DTAGrez(ideltaQLL13, deltaQLL13_)
	  DTAGimz(ideltaQLL13, deltaQLL13_)

	  DTAGrez(ideltaULR12, deltaULR12_)
	  DTAGimz(ideltaULR12, deltaULR12_)
	  DTAGrez(ideltaULR23, deltaULR23_)
	  DTAGimz(ideltaULR23, deltaULR23_)
	  DTAGrez(ideltaULR13, deltaULR13_)
	  DTAGimz(ideltaULR13, deltaULR13_)

	  DTAGrez(ideltaURL12, deltaURL12_)
	  DTAGimz(ideltaURL12, deltaURL12_)
	  DTAGrez(ideltaURL23, deltaURL23_)
	  DTAGimz(ideltaURL23, deltaURL23_)
	  DTAGrez(ideltaURL13, deltaURL13_)
	  DTAGimz(ideltaURL13, deltaURL13_)

	  DTAGrez(ideltaURR12, deltaURR12_)
	  DTAGimz(ideltaURR12, deltaURR12_)
	  DTAGrez(ideltaURR23, deltaURR23_)
	  DTAGimz(ideltaURR23, deltaURR23_)
	  DTAGrez(ideltaURR13, deltaURR13_)
	  DTAGimz(ideltaURR13, deltaURR13_)

	  DTAGrez(ideltaDLR12, deltaDLR12_)
	  DTAGimz(ideltaDLR12, deltaDLR12_)
	  DTAGrez(ideltaDLR23, deltaDLR23_)
	  DTAGimz(ideltaDLR23, deltaDLR23_)
	  DTAGrez(ideltaDLR13, deltaDLR13_)
	  DTAGimz(ideltaDLR13, deltaDLR13_)

	  DTAGrez(ideltaDRL12, deltaDRL12_)
	  DTAGimz(ideltaDRL12, deltaDRL12_)
	  DTAGrez(ideltaDRL23, deltaDRL23_)
	  DTAGimz(ideltaDRL23, deltaDRL23_)
	  DTAGrez(ideltaDRL13, deltaDRL13_)
	  DTAGimz(ideltaDRL13, deltaDRL13_)

	  DTAGrez(ideltaDRR12, deltaDRR12_)
	  DTAGimz(ideltaDRR12, deltaDRR12_)
	  DTAGrez(ideltaDRR23, deltaDRR23_)
	  DTAGimz(ideltaDRR23, deltaDRR23_)
	  DTAGrez(ideltaDRR13, deltaDRR13_)
	  DTAGimz(ideltaDRR13, deltaDRR13_)
	endif

	sf_valid = 0
	tl_valid = 0
	higgs_valid = 0
	coup_valid = 0

	deltaSf_LL(1,2,3) = deltaQLL12_
	deltaSf_LL(2,3,3) = deltaQLL23_
	deltaSf_LL(1,3,3) = deltaQLL13_

	deltaSf_LR(1,2,3) = deltaULR12_
	deltaSf_LR(2,3,3) = deltaULR23_
	deltaSf_LR(1,3,3) = deltaULR13_

	deltaSf_RL(1,2,3) = deltaURL12_
	deltaSf_RL(2,3,3) = deltaURL23_
	deltaSf_RL(1,3,3) = deltaURL13_

	deltaSf_RR(1,2,3) = deltaURR12_
	deltaSf_RR(2,3,3) = deltaURR23_
	deltaSf_RR(1,3,3) = deltaURR13_

	call IniFV(3)

	deltaSf_LL(1,2,4) = deltaQLL12_
	deltaSf_LL(2,3,4) = deltaQLL23_
	deltaSf_LL(1,3,4) = deltaQLL13_

	deltaSf_LR(1,2,4) = deltaDLR12_
	deltaSf_LR(2,3,4) = deltaDLR23_
	deltaSf_LR(1,3,4) = deltaDLR13_

	deltaSf_RL(1,2,4) = deltaDRL12_
	deltaSf_RL(2,3,4) = deltaDRL23_
	deltaSf_RL(1,3,4) = deltaDRL13_

	deltaSf_RR(1,2,4) = deltaDRR12_
	deltaSf_RR(2,3,4) = deltaDRR23_
	deltaSf_RR(1,3,4) = deltaDRR13_

	call IniFV(4)
	end


************************************************************************

	subroutine IniFV(t)
	implicit none
	integer t

#include "FH.h"

	integer nQ, nU
	RealType z

	nQ = tQ(t)
	MSS2(1,2,nQ) = deltaSf_LL(1,2,t)*MSS2_LL(1,2)
	z = Sq(MSS2(1,2,nQ))
	MSS2(2,3,nQ) = deltaSf_LL(2,3,t)*MSS2_LL(2,3)
	z = z + Sq(MSS2(2,3,nQ))
	MSS2(1,3,nQ) = deltaSf_LL(1,3,t)*MSS2_LL(1,3)
	z = z + Sq(MSS2(1,3,nQ))

	nU = tU(t)
	MSS2(1,2,nU) = deltaSf_RR(1,2,t)*MSS2_RR(1,2)
	z = z + Sq(MSS2(1,2,nU))
	MSS2(2,3,nU) = deltaSf_RR(2,3,t)*MSS2_RR(2,3)
	z = z + Sq(MSS2(2,3,nU))
	MSS2(1,3,nU) = deltaSf_RR(1,3,t)*MSS2_RR(1,3)
	z = z + Sq(MSS2(1,3,nU))

	Kf(1,2,t) = deltaSf_LR(1,2,t)*MSS2_LR(1,2)
	z = z + Sq(Kf(1,2,t))
	Kf(2,3,t) = deltaSf_LR(2,3,t)*MSS2_LR(2,3)
	z = z + Sq(Kf(2,3,t))
	Kf(1,3,t) = deltaSf_LR(1,3,t)*MSS2_LR(1,3)
	z = z + Sq(Kf(1,3,t))
	Kf(2,1,t) = Conjugate(deltaSf_RL(1,2,t)*MSS2_RL(1,2))
	z = z + Sq(Kf(2,1,t))
	Kf(3,2,t) = Conjugate(deltaSf_RL(2,3,t)*MSS2_RL(2,3))
	z = z + Sq(Kf(3,2,t))
	Kf(3,1,t) = Conjugate(deltaSf_RL(1,3,t)*MSS2_RL(1,3))
	z = z + Sq(Kf(3,1,t))

	if( z .lt. 1D-13 ) then
	  fv = ibclr(fv, t)
	else
	  fv = ibset(fv, t)
	  if( debuglevel .ge. 1 )
     &      DPARA "NMFV on for Sfermion ", Digit(t) ENDL
	endif
	end

