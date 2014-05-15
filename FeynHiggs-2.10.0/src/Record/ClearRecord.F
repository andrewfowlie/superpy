* ClearRecord.F
* initialize a FH record
* this file is part of FeynHiggs
* last modified 9 Aug 11 th

#include "externals.h"
#include "types.h"


	subroutine FHClearRecord(record)
	implicit none

#include "FHRecord.h"

	RecordDecl(record)

	integer v

	do v = 1, FHRecordC - 1
	  record(v,iVar) = 0
	  record(v,iLower) = unset
	  record(v,iUpper) = 0
	  record(v,iStep) = 0
	enddo

	do v = FHRangeC
	  record(iRe(v),iVar) = 0
	  record(iRe(v),iLower) = unset
	  record(iRe(v),iUpper) = 0
	  record(iRe(v),iStep) = 0

	  record(iIm(v),iVar) = 0
	  record(iIm(v),iLower) = 0
	  record(iIm(v),iUpper) = 0
	  record(iIm(v),iStep) = 0

	  record(iAbs(v),iVar) = 0
	  record(iAbs(v),iLower) = unset
	  record(iAbs(v),iUpper) = 0
	  record(iAbs(v),iStep) = 0

	  record(iArg(v),iVar) = 0
	  record(iArg(v),iLower) = 0
	  record(iArg(v),iUpper) = 0
	  record(iArg(v),iStep) = 0
	enddo
	end

