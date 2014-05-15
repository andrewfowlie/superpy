* LoopRecord.F
* iterate the FH record
* this file is part of FeynHiggs
* last modified 30 Nov 11 th

#include "externals.h"
#include "types.h"
#include "debug.h"

#define SKIP "% "
#define Info print *, SKIP,

#define VERBOSE record(iAdmin,iStep) .ne. 0

* the fields record(iAdmin,*) contain the following
* administrative entries (integer):
*   record(iAdmin,iVar)   = loop count (0 if not yet initialized),
*   record(iAdmin,iLower) = input field #1 for interpolating a table
*   record(iAdmin,iUpper) = input field #2 for interpolating a table
*   record(iAdmin,iStep)  = verbosity of FHLoopRecord

	subroutine FHLoopRecord(error, record)
	implicit none
	integer error

#include "FHRecord.h"
#include "const.h"
#define __SUBROUTINE__ "FHLoopRecord"

	RecordDecl(record)

	integer v
	ComplexType tmpc
	RealType tmpr

#if U77EXT
#define len_trim lnblnk
	integer lnblnk
	external lnblnk
#endif

	error = 0

	if( record(iAdmin,iVar) .eq. 0 ) then
* FIRST ITERATION

* do initial assignment of defaults

#define RConstDefault(v,val) call RIniDefault(record,v,0,val)
#define RVarDefault(v,d) call RIniDefault(record,v,d,tmpr)
#define CConstDefault(v,val) call CIniDefault(record,v,0,val)
#define CVarDefault(v,d) call CIniDefault(record,v,d,tmpc)

#include "Defaults.h"

#undef RConstDefault
#undef CConstDefault
#undef RVarDefault
#undef CVarDefault

* check for valid arguments

	  do v = 0, 8, 4
	    if( record(iRe(iXtau+v),iLower) .ne. unset .or.
     &          record(iAbs(iXtau+v),iLower) .ne. unset ) then
	      if( VERBOSE .and.
     &            (record(iRe(iAtau+v),iLower) .ne. unset .or.
     &             record(iAbs(iAtau+v),iLower) .ne. unset) )
     &          Info Strip(FHName(iXtau+v)),
     &            " supersedes ", Strip(FHName(iAtau+v))
	      record(iRe(iAtau+v),iLower) = -1
	    else
	      record(iRe(iXtau+v),iLower) = -1
	    endif
	  enddo

	  if( record(iMA0,iLower) .ge. 0 ) then
	    if( VERBOSE .and. record(iMHp,iLower) .ge. 0 )
     &        Info "MA0 supersedes MHp"
	    record(iMHp,iLower) = -1
	  else
	    record(iMA0,iLower) = -1
	  endif

	  do v = FHRangeR
	    if( record(v,iLower) .eq. unset .and.
     &        (v .ne. iMSusy .or. record(iM3SQ,iLower) .eq. default) )
     &        Error(error, Strip(FHName(v))//" not given")
	  enddo

	  do v = FHRangeC
	    if( record(iRe(v),iLower) .eq. unset .and.
     &          record(iAbs(v),iLower) .eq. unset )
     &        Error(error, "Neither Re nor Abs("//FHNameC(v)//") given")
	  enddo

	  do v = FHRangeA
	    record(v,iVar) = record(v,iLower)
	  enddo

	else
* NEXT ITERATION

	  do v = FHRangeA
	    if( record(v,iStep) .ne. 0 ) then
	      if( record(v,iStep)*(record(v,iUpper) -
     &                             record(v,iLower)) .lt. 0 ) then
	        record(v,iVar) = record(v,iVar)*abs(record(v,iStep))
	      else
	        record(v,iVar) = record(v,iVar) + record(v,iStep)
	      endif
	      if( record(v,iVar) - record(v,iUpper) .lt. 1D-5 ) goto 1
	      record(v,iVar) = record(v,iLower)
	    endif
	  enddo

	  error = -1
	  return

1	  continue

	endif

	if( record(iAdmin,iLower) .ne. unset ) then
	  call Interpolate(error, record)
	  if( error .ne. 0 ) return
	endif

	do v = FHRangeC
	  if( record(iAbs(v),iVar) .ne. unset ) then
	    tmpc = record(iAbs(v),iVar)*
     &        exp(cI*record(iArg(v),iVar))
	    record(iRe(v),iVar) = Re(tmpc)
	    record(iIm(v),iVar) = Im(tmpc)
	  endif
	enddo

#define MUExTB(h,op) record(h(iMUE),iVar) op record(iTB,iVar)

	if( record(iRe(iXtau),iVar) .ne. -1 ) then
	  record(iRe(iAtau),iVar) = record(iRe(iXtau),iVar) + MUExTB(iRe, *)
	  record(iIm(iAtau),iVar) = record(iIm(iXtau),iVar) - MUExTB(iIm, *)
	endif
	if( record(iRe(iXt),iVar) .ne. -1 ) then
	  record(iRe(iAt),iVar) = record(iRe(iXt),iVar) + MUExTB(iRe, /)
	  record(iIm(iAt),iVar) = record(iIm(iXt),iVar) - MUExTB(iIm, /)
	endif
	if( record(iRe(iXb),iVar) .ne. -1 ) then
	  record(iRe(iAb),iVar) = record(iRe(iXb),iVar) + MUExTB(iRe, *)
	  record(iIm(iAb),iVar) = record(iIm(iXb),iVar) - MUExTB(iIm, *)
	endif

* do final assignment of `inherited' values

#define RConstDefault(v,val)
#define CConstDefault(v,val)
#define RVarDefault(v,d) call RSetDefault(record, v, d)
#define CVarDefault(v,d) call CSetDefault(record, v, d)

#include "Defaults.h"

#undef RConstDefault
#undef CConstDefault
#undef RVarDefault
#undef CVarDefault

	record(iAdmin,iVar) = record(iAdmin,iVar) + 1
	end


************************************************************************

	subroutine RIniDefault(record, v, d, val)
	implicit none
	integer v, d
	RealType val

#include "FHRecord.h"

	RecordDecl(record)

#if U77EXT
	integer lnblnk
	external lnblnk
#endif

	integer i, j

	if( record(v,iLower) .ne. unset ) return

	if( d .eq. 0 ) then
	  if( VERBOSE )
     &      Info FHNameR(v), " defaults to ", val
	  record(v,iLower) = val
	else
	  if( VERBOSE ) then
	    i = d
1	    j = iand(i, 255)
	    i = ishft(i, -8)
	    if( i .ne. 0 .and. record(j,iLower) .eq. unset ) goto 1
	    Info FHNameR(v), " defaults to ", FHNameR(j)
	  endif
	  record(v,iLower) = default
	endif
	end


************************************************************************

	subroutine CIniDefault(record, v, d, val)
	implicit none
	integer v, d
	ComplexType val

#include "FHRecord.h"

	RecordDecl(record)

#if U77EXT
	integer lnblnk
	external lnblnk
#endif

	integer i, j

	if( record(iAbs(v),iLower) .ne. unset ) then
	  if( VERBOSE .and. record(iRe(v),iLower) .ne. unset )
     &      Info FHNameC(v), ": Abs/Arg supersedes Re/Im"
	  return
	endif

	if( record(iRe(v),iLower) .ne. unset ) return

	if( d .eq. 0 ) then
	  if( VERBOSE )
     &      Info FHNameC(v), " defaults to ", val
	  record(iRe(v),iLower) = Re(val)
	  record(iIm(v),iLower) = Im(val)
	else
	  if( VERBOSE ) then
	    i = d
1	    j = iand(i, 255)
	    i = ishft(i, -8)
	    if( i .ne. 0 .and. record(iRe(j),iLower) .eq. unset ) goto 1
	    Info FHNameC(v), " defaults to ", FHNameC(j)
	  endif
	  record(iRe(v),iLower) = default
	endif
	end


************************************************************************

	subroutine RSetDefault(record, v, d)
	implicit none
	integer v, d

#include "FHRecord.h"

	RecordDecl(record)

	integer i, j

	i = d
	do while( i .ne. 0 )
	  j = iand(i, 255)
	  if( record(v,iLower) .eq. default .and.
     &        record(j,iLower) .ne. unset ) then
	    record(v,iVar) = record(j,iVar)
	    return
	  endif
	  i = ishft(i, -8)
	enddo
	end


************************************************************************

	subroutine CSetDefault(record, v, d)
	implicit none
	integer v, d

#include "FHRecord.h"

	RecordDecl(record)

	integer i, j

	i = d
	do while( i .ne. 0 )
	  j = iand(i, 255)
	  if( record(iRe(v),iLower) .eq. default .and.
     &        (record(iRe(j),iLower) .ne. unset .or.
     &         record(iAbs(j),iLower) .ne. unset) ) then
	    record(iRe(v),iVar) = record(iRe(j),iVar)
	    record(iIm(v),iVar) = record(iIm(j),iVar)
	    return
	  endif
	  i = ishft(i, -8)
	enddo
	end


************************************************************************

	subroutine Interpolate(error, record)
	implicit none
	integer error

#include "FH.h"
#include "FHRecord.h"

	RecordDecl(record)

	integer nfit
	RealType eps
	parameter (nfit = 4)
	parameter (eps = 1D-6)

	integer i1, i2, v, f, r, o, missing
	RealType val1, val2, cur1, cur2
	RealType ref1, ref2, d1, d2, a1, a2
	RealType vdm(nfit,nfit), b(nfit)
	integer row(nfit), perm(nfit)

	integer compass
	RealType e1, e2
	compass(e1, e2) = int(sign(.5D0, d1 - e1) +
     &    sign(1D0, d2 - e2) + 2.5D0)

	error = 0

	i1 = int(record(iAdmin,iLower))
	i2 = int(record(iAdmin,iUpper))
	val1 = record(i1,iVar)
	val2 = record(i2,iVar)

	ref1 = 0
	ref2 = 0
	a1 = 1D300
	a2 = 1D300

	do f = 1, nfit
	  row(f) = 0
	enddo
	missing = 30

	do r = 1, tablerows
	  cur1 = tabledata(i1,r)
	  cur2 = tabledata(i2,r)

	  d1 = abs(cur1 - ref1)
	  if( d1 .gt. 1D-13 .and. d1 .lt. a1 ) then
	    a1 = d1
	    ref1 = cur1
	  endif
	  d2 = abs(cur2 - ref2)
	  if( d2 .gt. 1D-13 .and. d2 .lt. a2 ) then
	    a2 = d2
	    ref2 = cur2
	  endif

	  d1 = cur1 - val1
	  d2 = cur2 - val2
	  f = compass(0D0, 0D0)
	  o = row(f)
	  if( o .eq. 0 .or.
     &       (abs(d1) .le. abs(tabledata(i1,o) - val1) .and.
     &        abs(d2) .le. abs(tabledata(i2,o) - val2)) ) then
	    row(f) = r
	    missing = ibclr(missing, f)
	  endif
	enddo

* catch points exactly on the border
	if( missing .ne. 0 ) then
	  do f = 1, nfit
	    r = row(f)
	    if( r .ne. 0 ) then
	      d1 = tabledata(i1,r) - val1
	      d2 = tabledata(i2,r) - val2
	      row(compass(eps, 0D0)) = r
	      row(compass(0D0, eps)) = r
	      row(compass(eps, eps)) = r
	    endif
	  enddo
	endif

	if( debuglevel .ge. 2 ) then
	  DPARA "Interpolation points chosen:" ENDL
	  DPARA Letter(1), tabledata(i1,row(1)), tabledata(i2,row(1)) ENDL
	  DPARA Letter(2), tabledata(i1,row(2)), tabledata(i2,row(2)) ENDL
	  DPARA Letter(3), tabledata(i1,row(3)), tabledata(i2,row(3)) ENDL
	  DPARA Letter(4), tabledata(i1,row(4)), tabledata(i2,row(4)) ENDL
	  DPARA "Grid spacing:" ENDL
	  DPARA a1, a2 ENDL
	endif

	a1 = 1.3D0*a1
	a2 = 1.3D0*a2

	do f = 1, nfit
	  r = row(f)
	  cur1 = tabledata(i1,r)
	  cur2 = tabledata(i2,r)

	  if( cur1 .eq. unset .or.
     &        cur2 .eq. unset .or.
     &        abs(cur1 - val1) .gt. a1 .or.
     &        abs(cur2 - val2) .gt. a2 ) then
	    Error(error, "Extrapolation not possible")
	    error = -2
	    return
	  endif

	  vdm(1,f) = 1
	  vdm(2,f) = cur1
	  vdm(3,f) = cur2
	  vdm(4,f) = cur1*cur2
	enddo

	call LUDecompose(nfit, vdm,nfit, perm)

	do v = FHRangeA
	  if( v .ne. i1 .and. v .ne. i2 .and.
     &        tableflag(v) .ne. 0 ) then
	    do f = 1, nfit
	      b(perm(f)) = tabledata(v,row(f))
	      if( debuglevel .ge. 4 )
     &          DPARA FHName(v), Letter(f), " ", b(perm(f)) ENDL
	    enddo
	    call LUSolve(nfit, vdm,nfit, b)
	    record(v,iVar) = b(1) + val1*b(2) + val2*(b(3) + val1*b(4))
	    if( debuglevel .ge. 3 )
     &        DPARA FHName(v), record(v,iVar) ENDL
	  endif
	enddo
	end

