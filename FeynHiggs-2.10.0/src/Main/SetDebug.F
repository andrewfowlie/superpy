* SetDebug.F
* set the debug level
* this file is part of FeynHiggs
* last modified 9 Aug 11 th

#include "externals.h"
#include "types.h"


	subroutine FHSetDebug(debuglevel_)
	implicit none
	integer debuglevel_

#include "FH.h"

	debuglevel = debuglevel_
	end

