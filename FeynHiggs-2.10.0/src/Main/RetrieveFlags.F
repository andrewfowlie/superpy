* RetrieveFlags.F
* retrieve all flags from FeynHiggs
* this file is part of FeynHiggs
* last modified 9 Aug 11 th

#include "externals.h"
#include "types.h"
#include "debug.h"


	subroutine FHRetrieveFlags(error,
     &    mssmpart_, fieldren_, tanbren_,
     &    higgsmix_, p2approx_, looplevel_,
     &    runningMT_, botResum_, tlCplxApprox_)
	implicit none
	integer error
	integer mssmpart_, fieldren_, tanbren_
	integer higgsmix_, p2approx_, looplevel_
	integer runningMT_, botResum_, tlCplxApprox_

#include "FH.h"
#define __SUBROUTINE__ "FHRetrieveFlags"

	error = 0
	CheckFlags(error)

	mssmpart_ = mssmpart
	fieldren_ = fieldren
	tanbren_ = tanbren
	higgsmix_ = higgsmix
	p2approx_ = p2approx
	looplevel_ = looplevel
	runningMT_ = runningMT
	botResum_ = botResum
	tlCplxApprox_ = tlCplxApprox
	end

#undef __SUBROUTINE__


************************************************************************

	subroutine FHRetrieveFlagsString(error, flags)
	implicit none
	integer error
	character*(*) flags

#include "FH.h"
#define __SUBROUTINE__ "FHRetrieveFlagsString"

	error = 0

	if( len(flags) .lt. 9 ) then
	  Error(error, "string not wide enough")
	  return
	endif

	CheckFlags(error)

	flags(1:1) = Digit(mssmpart)
	flags(2:2) = Digit(fieldren)
	flags(3:3) = Digit(tanbren)
	flags(4:4) = Digit(higgsmix)
	flags(5:5) = Digit(p2approx)
	flags(6:6) = Digit(looplevel)
	flags(7:7) = Digit(runningMT)
	flags(8:8) = Digit(botResum)
	flags(9:9) = Digit(tlCplxApprox)
        end

