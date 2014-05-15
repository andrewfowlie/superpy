* RecordIndex.F
* convert variable names to indices
* this file is part of FeynHiggs
* last modified 15 Apr 12 th

#include "externals.h"
#include "types.h"
#include "debug.h"


	subroutine FHRecordIndex(ind, para)
	implicit none
	integer ind
	character*(*) para

#include "FHRecord.h"
#define __SUBROUTINE__ "FHRecordIndex"

	integer v
	character*128 line

#if U77EXT
#define len_trim lnblnk
	integer lnblnk
	external lnblnk
#endif

	external recordnamesini

	do v = FHRangeA
	  if( para .eq. FHName(v) ) then
	    ind = v
	    return
	  endif
	enddo

* e.g. At is shorthand for Re(At)
	do v = FHRangeC
	  if( para .eq. FHNameC(v) ) then
	    ind = v
	    return
	  endif
	enddo

	line = "Unknown parameter "//para
	Warning(Strip(line))
	ind = 0
	end

