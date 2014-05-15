* TableRecord.F
* associate a table with an FH record
* this file is part of FeynHiggs
* last modified 9 Aug 11 th

#include "externals.h"
#include "types.h"
#include "debug.h"


	subroutine FHTableRecord(error, record, i1, i2)
	implicit none
	integer error, i1, i2

#include "FH.h"
#include "FHRecord.h"
#define __SUBROUTINE__ "FHTableRecord"

	RecordDecl(record)

	integer v

	error = 0

	if( i1 .lt. FHRecordR .or. i1 .ge. FHRecordE .or.
     &      tableflag(i1) .eq. 0 .or.
     &      i2 .lt. FHRecordR .or. i2 .ge. FHRecordE .or.
     &      tableflag(i2) .eq. 0 ) then
	  Error(error, "Invalid interpolation point")
	  return
	endif

	record(iAdmin,iLower) = i1
	record(iAdmin,iUpper) = i2

	do v = FHRangeA
	  if( v .ne. i1 .and. v .ne. i2 .and.
     &        tableflag(v) .ne. 0 )
     &      record(v,iLower) = bytable
	enddo
	end

