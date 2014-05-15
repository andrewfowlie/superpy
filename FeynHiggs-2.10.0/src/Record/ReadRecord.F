* ReadRecord.F
* read a parameter file into a FH record
* this file is part of FeynHiggs
* last modified 30 Nov 11 th

#include "externals.h"
#include "types.h"
#include "debug.h"


	subroutine FHReadRecord(error, record, slhadata, inputfile)
	implicit none
	integer error
	character*(*) inputfile

#include "FH.h"
#include "FHRecord.h"
#include "SLHA.h"
#define __SUBROUTINE__ "FHReadRecord"

	RecordDecl(record)
	ComplexType slhadata(nslhadata)

	RealType step
	integer unit, token, last, v, ind1, ind2
	logical islog, isdiv

	integer NextToken
	external NextToken

#if U77EXT
#define len_trim lnblnk
	integer lnblnk
	external lnblnk
#endif

	integer lineno
	character*200 line
	common /slhaline/ lineno, line

	character*(*) incompat
	parameter (incompat =
     &    "Opposite-sign limits incompatible with log step for ")

	error = 0

	call SLHARead(error, slhadata, inputfile, 1)
	if( error .eq. 0 ) then
	  call FHSLHARecord(error, record, slhadata)
	  return
	endif

	call FHClearRecord(record)

* continue if only not an SLHA file
	if( error .ne. 2 ) return

	if( inputfile .eq. "-" ) then
	  unit = 5
	  if( lineno .ne. 0 ) goto 101
	else
	  open(10, file=inputfile, status="old", err=97)
	  unit = 10
	endif

	lineno = 0

1	continue
	lineno = lineno + 1

100	format(A)
	read(unit, 100, err=99, end=99) line

101	last = 1

	token = NextToken(line, last)
	if( token .eq. 0 ) goto 1

	if( line(token:last-1) .eq. "table" ) then
	  token = NextToken(line, last)
	  if( token .eq. 0 ) goto 98
	  call FHLoadTable(error, line(token:last-1), unit)
	  if( error .ne. 0 ) return

	  token = NextToken(line, last)
	  if( token .eq. 0 ) goto 98
	  call FHRecordIndex(ind1, line(token:last-1))

	  token = NextToken(line, last)
	  if( token .eq. 0 ) goto 98
	  call FHRecordIndex(ind2, line(token:last-1))

	  call FHTableRecord(error, record, ind1, ind2)
	  if( error .ne. 0 ) return
	  goto 1
	endif

	call FHRecordIndex(v, line(token:last-1))
	if( v .eq. 0 ) goto 1

	record(v,iLower) = unset
	record(v,iStep) = 0

	token = NextToken(line, last)
	if( token .eq. 0 ) goto 98
	read(line(token:last-1), *, err=98, end=98) record(v,iLower)
	if( debuglevel .ge. 2 )
     &    DPARA FHName(v), " lower =", record(v,iLower) ENDL

	token = NextToken(line, last)
	if( token .eq. 0 ) goto 1
	read(line(token:last-1), *, err=98, end=98) record(v,iUpper)
	record(v,iStep) = sign(1D0, record(v,iUpper) - record(v,iLower))
	if( debuglevel .ge. 2 )
     &    DPARA FHName(v), " upper =", record(v,iUpper) ENDL

	token = NextToken(line, last)
	if( token .eq. 0 ) goto 1

	islog = line(token:token) .eq. '*'
	if( islog ) then
	  last = token + 1
	  token = NextToken(line, last)
	  if( token .eq. 0 ) goto 98
	endif

	isdiv = line(token:token) .eq. '/'
	if( isdiv ) then
	  last = token + 1
	  token = NextToken(line, last)
	  if( token .eq. 0 ) goto 98
	endif
	read(line(token:last-1), *, err=98, end=98) step
	if( debuglevel .ge. 2 )
     &    DPARA FHName(v), "  step =", step ENDL

	if( islog .and. record(v,iLower)*record(v,iUpper) .le. 0 ) then
	  Error(error, incompat//Strip(FHName(v)))
	  return
	endif

	if( isdiv ) then
	  step = 1/max(step - 1, 1D0)
	  if( islog ) then
	    step = (record(v,iUpper)/record(v,iLower))**step
	  else
	    step = (record(v,iUpper) - record(v,iLower))*step
	  endif
	endif

	if( islog ) then
	  record(v,iStep) = -sign(step, record(v,iStep))
	else
	  if( step*record(v,iStep) .lt. 0 ) step = 0
	  record(v,iStep) = step
	endif

	if( NextToken(line, last) .eq. 0 ) goto 1

998	format("Syntax error in ", A, ", line", I3)
98	write(line, 998) Strip(inputfile), lineno
	Warning(Strip(line))
	goto 1

97	line = Strip(inputfile)//" not found"
	Error(error, Strip(line))
	return

99	close(unit)
* indicate non-SLHA file (same code as SLHARead)
	error = 2
	end

