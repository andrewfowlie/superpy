* Flavour.F
* evaluation of flavour observables
* used as further constraints for the MSSM parameter space
* this file is part of FeynHiggs
* last modified 13 Jun 12 th

#include "externals.h"
#include "types.h"
#include "debug.h"


	subroutine FHFlavour(error,
     &    bsgMSSM, bsgSM,
     &    deltaMsMSSM, deltaMsSM,
     &    bsmumuMSSM, bsmumuSM)
	implicit none
	integer error
	RealType bsgMSSM, bsgSM
	RealType deltaMsMSSM, deltaMsSM
	RealType bsmumuMSSM, bsmumuSM

#include "gm2.h"
#define __SUBROUTINE__ "FHFlavour"

	call CheckHiggsCorr(error)
	if( error .gt. 0 ) return

* a) b -> s gamma
	call Bsgamma(bsgMSSM, bsgSM)

* b) delta Ms
	call DeltaMs(deltaMsMSSM, deltaMsSM)

* c) Bs -> mu+ mu-
	call Bsll(bsmumuMSSM, bsmumuSM)
c	bsmumuMSSM = -999
c	bsmumuSM = -999
	end

