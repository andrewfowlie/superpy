* ini.F
* routines for initializing and setting some parameters
* this file is part of FeynHiggs 
* last modified 30 Nov 11 th

#include "externals.h"
#include "types.h"
#include "ltdefs.h"


	subroutine clearcache
	implicit none

#include "lt.h"
#include "ltcache.h"

	integer*8 cacheptr(4,ncaches)
	equivalence (cache, cacheptr)

	integer i

	do i = 1, ncaches
	  cacheptr(1,i) = 0
	  cacheptr(2,i) = 0
	enddo
	end

************************************************************************

	subroutine setmudim(mudim_)
	implicit none
	RealType mudim_

#include "lt.h"

	mudim = mudim_
	call clearcache
	end

************************************************************************

	RealType function getmudim()
	implicit none

#include "lt.h"

	getmudim = mudim
	end

************************************************************************

	subroutine setdelta(delta_)
	implicit none
	RealType delta_

#include "lt.h"

	delta = delta_
	call clearcache
	end

************************************************************************

	RealType function getdelta()
	implicit none

#include "lt.h"

	getdelta = delta
	end

************************************************************************

	subroutine setlambda(lambda_)
	implicit none
	RealType lambda_

#include "lt.h"

	lambda = lambda_
	call clearcache
	end

************************************************************************

	RealType function getlambda()
	implicit none

#include "lt.h"

	getlambda = lambda
	end

************************************************************************

	subroutine ltini()
	implicit none

#include "lt.h"

	delta = 0
	mudim = 1
	lambda = 1
	call clearcache
	end

************************************************************************

	subroutine ltexi()
	end

************************************************************************

	block data ltcutoffini

#include "lt.h"

	data delta /0/
	data mudim /1/
	data lambda /1/
	end

