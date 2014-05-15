* RetrieveSMPara.F
* retrieve the SM parameters
* this file is part of FeynHiggs
* last modified 30 Nov 11 th

#include "externals.h"
#include "types.h"
#include "debug.h"


	subroutine FHRetrieveSMPara(error,
     &    invAlfaMZ_, AlfasMZ_, GF_,
     &    ME_, MU_, MD_, MM_, MC_, MS_, ML_, MB_,
     &    MW_, MZ_,
     &    CKMlambda_, CKMA_, CKMrhobar_, CKMetabar_)
	implicit none
	integer error
	RealType invAlfaMZ_, AlfasMZ_, GF_
	RealType ME_, MU_, MD_, MM_, MC_, MS_, ML_, MB_
	RealType MW_, MZ_
	RealType CKMlambda_, CKMA_, CKMrhobar_, CKMetabar_

#include "FH.h"
#define __SUBROUTINE__ "FHRetrieveSMPara"

	error = 0
	CheckSMPara(error)

	invAlfaMZ_ = invAlfaMZ
	AlfasMZ_ = AlfasMZ
	GF_ = GF

	ME_ = ME
	MU_ = MU
	MD_ = MD
	MM_ = MM
	MC_ = MC
	MS_ = MS
	ML_ = ML
	MB_ = MB

	MW_ = MW
	MZ_ = MZ

	CKMlambda_ = CKMlambda
	CKMA_ = CKMA
	CKMrhobar_ = CKMrhobar
	CKMetabar_ = CKMetabar
	end

