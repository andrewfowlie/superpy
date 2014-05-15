* SetSLHA.F
* initialize parameters from the SLHA data structure
* this file is part of FeynHiggs
* last modified 30 Nov 11 th

#include "externals.h"
#include "types.h"
#include "debug.h"


	subroutine FHSetSLHA(error, slhadata)
	implicit none
	integer error
	ComplexType slhadata(*)

#include "SLHA.h"
#include "FH.h"
#define __SUBROUTINE__ "FHSetSLHA"

#include "FHRecord.h"

	RecordDecl(record)

	call FHSLHARecord(error, record, slhadata)
	if( error .ne. 0 ) return

	call FHLoopRecord(error, record)
	if( error .ne. 0 ) return

	call FHSetRecord(error, record)
	end

