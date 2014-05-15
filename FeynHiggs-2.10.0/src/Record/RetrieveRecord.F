* RetrieveRecord.F
* initialize the record from currently set FH parameters
* this file is part of FeynHiggs
* last modified 29 May 12 th

#include "externals.h"
#include "types.h"
#include "debug.h"


	subroutine FHRetrieveRecord(error, record, iX)
	implicit none
	integer error, iX

#include "FHRecord.h"
#define __SUBROUTINE__ "FHRetrieveRecord"

	RecordDecl(record)

#define R(i) record(i,iX)

	call FHRetrieveSMPara(error,
     &    R(iinvAlfaMZ), R(iAlfasMZ), R(iGF),
     &    R(iME), R(iMU), R(iMD),
     &    R(iMM), R(iMC), R(iMS),
     &    R(iML),         R(iMB),
     &    R(iMW), R(iMZ),
     &    R(iCKMlambda), R(iCKMA),
     &    R(iCKMrhobar), R(iCKMetabar))
	if( error .ne. 0 ) return

	call FHRetrievePara(error, R(iscalefactor), R(iMT),
     &    R(iTB), R(iMA0), R(iMHp),
     &    R(iM3SL), R(iM3SE),
     &    R(iM3SQ), R(iM3SU), R(iM3SD),
     &    R(iM2SL), R(iM2SE),
     &    R(iM2SQ), R(iM2SU), R(iM2SD),
     &    R(iM1SL), R(iM1SE),
     &    R(iM1SQ), R(iM1SU), R(iM1SD),
     &    R(iMUE),
     &    R(iAtau), R(iAt), R(iAb),
     &    R(iAmu), R(iAc), R(iAs),
     &    R(iAe), R(iAu), R(iAd),
     &    R(iM1), R(iM2), R(iM3),
     &    R(iQtau), R(iQt), R(iQb))
	if( error .ne. 0 ) return

	call FHRetrieveNMFV(error,
     &    R(ideltaQLL12), R(ideltaQLL23), R(ideltaQLL13),
     &    R(ideltaULR12), R(ideltaULR23), R(ideltaULR13),
     &    R(ideltaURL12), R(ideltaURL23), R(ideltaURL13),
     &    R(ideltaURR12), R(ideltaURR23), R(ideltaURR13),
     &    R(ideltaDLR12), R(ideltaDLR23), R(ideltaDLR13),
     &    R(ideltaDRL12), R(ideltaDRL23), R(ideltaDRL13),
     &    R(ideltaDRR12), R(ideltaDRR23), R(ideltaDRR13))

	call FHRetrieveLFV(error,
     &    R(ideltaLLL12), R(ideltaLLL23), R(ideltaLLL13),
     &    R(ideltaELR12), R(ideltaELR23), R(ideltaELR13),
     &    R(ideltaERL12), R(ideltaERL23), R(ideltaERL13),
     &    R(ideltaERR12), R(ideltaERR23), R(ideltaERR13))
	end

