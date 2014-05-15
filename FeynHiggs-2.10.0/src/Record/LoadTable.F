* LoadTable.F
* load a table file into memory
* this file is part of FeynHiggs
* last modified 9 Aug 11 th

#include "externals.h"
#include "types.h"
#include "debug.h"


	subroutine FHLoadTable(error, inputfile, inputunit)
	implicit none
	integer error, inputunit
	character*(*) inputfile

#include "FHRecord.h"
#define __SUBROUTINE__ "FHLoadTable"

	character*1024 line
	integer unit, col, ncol, token, last, ios
	integer ind(2*maxcols)

	integer NextToken
	external NextToken

#if U77EXT
#define len_trim lnblnk
	integer lnblnk
	external lnblnk
#endif

	error = 0
	tablerows = 0

	if( inputfile .eq. "-" ) then
	  unit = inputunit
	else
	  open(11, file=inputfile, status="old", err=97)
	  unit = 11
	endif

	do col = 1, maxcols
	  tableflag(col) = 0
	enddo

100	format(A)
	read(unit, 100, end=99) line

	last = 1
	do ncol = 1, 2*maxcols
	  token = NextToken(line, last)
	  if( token .eq. 0 ) goto 1
	  call FHRecordIndex(ind(ncol), line(token:last-1))
	  tableflag(ind(ncol)) = 1
	enddo

1	do tablerows = 0, maxrows - 1
	  read(unit, 100, end=99) line

	  last = 1
	  do col = 1, ncol - 1
	    token = NextToken(line, last)
	    if( token .eq. 0 ) goto 2

	    if( ind(col) .ne. 0 ) then
	      read(line(token:last-1), *, iostat=ios)
     &          tabledata(ind(col),tablerows+1)
	      if( ios .ne. 0 ) then
101	        format("Syntax error in ", A, ", line", I3)
	        write(line, 101) Strip(inputfile), tablerows + 2
	        Warning(Strip(line))
	      endif
	    endif
	  enddo

2	  continue
	enddo

99	close(unit)
	return

97	line = Strip(inputfile)//" not found"
	Error(error, Strip(line))
	end

