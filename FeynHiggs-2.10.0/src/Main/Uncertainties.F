* Uncertainties.F
* Compute (actually, not completely invent) an error estimate
* for the Higgs masses and mixings
* this file is part of FeynHiggs
* last modified 10 May 12 th

#include "externals.h"
#include "types.h"
#include "debug.h"


	subroutine FHUncertainties(error,
     &    DeltaMHiggs_, DeltaSAeff_, DeltaUHiggs_, DeltaZHiggs_)
	implicit none

#include "FH.h"
#include "looptools.h"
#define __SUBROUTINE__ "FHUncertainties"

	integer error
	RealType DeltaMHiggs_(NHiggs)
	ComplexType DeltaSAeff_
	ComplexType DeltaUHiggs_(NNeutral,NNeutral)
	ComplexType DeltaZHiggs_(NNeutral,NNeutral)

	RealType HiggsErr(NHiggsErr)
	RealType HiggsSave(NHiggsData)
	RealType HiggsBase(NHiggsErr)
	RealType oldmudim
	integer i, j, err

	integer NSfData
	parameter (NSfData = NSf - 9)
* Sf(*,bTR0) minus the CKM copy in bTR0

	RealType SfSave(NSfData)

	RealType SfData(NSfData)
	equivalence (Sf(1,bTR0), SfData)

	RealType DeltaMHiggs(NHiggs)
	ComplexType DeltaSAeff
	ComplexType DeltaXHiggs(0:NNeutral,0:NNeutral,2)
	equivalence (HiggsErr(1), DeltaMHiggs)
	equivalence (HiggsErr(1+NHiggs), DeltaSAeff)
	equivalence (HiggsErr(1+NHiggs+2), DeltaXHiggs)

	do i = 1, NHiggs
	  DeltaMHiggs_(i) = 100000
	enddo
	DeltaSAeff_ = uncomputable
	do i = 1, NNeutral
	  do j = 1, NNeutral
	    DeltaUHiggs_(j,i) = uncomputable
	    DeltaZHiggs_(j,i) = uncomputable
	  enddo
	enddo

	if( looplevel .lt. 2 ) return

	call CheckHiggsCorr(error)

	if( error .gt. 0 ) return

* store the central value

	do i = 1, NHiggsData
	  HiggsSave(i) = HiggsData(i)
	enddo

	do i = 1, NHiggsErr
	  HiggsBase(i) = HiggsData(i)
	  HiggsErr(i) = uncomputable
	enddo

	err = 0

* a) variation of mudim from MT/2 to 2 MT

	if( debuglevel .ge. 2 )
     &    DHIGGS "Uncertainties: mudim = 1/2 MT" ENDL

	oldmudim = getmudim()

	call setmudim(.25D0*MT2)
	call CalcHiggsCorr(err)

	if( err .ne. 0 ) then
	  Warning("Could not evaluate corrections at mudim = MT/2")
	  call setmudim(oldmudim)
	  goto 999
	endif

	do i = 1, NHiggsErr
	  if( HiggsData(i) .ne. uncomputable .and.
     &        HiggsBase(i) .ne. uncomputable )
     &      HiggsErr(i) = abs(HiggsData(i) - HiggsBase(i))
	enddo

	if( debuglevel .ge. 2 )
     &    DHIGGS "Uncertainties: mudim = 2 MT" ENDL

	call setmudim(4*MT2)
	call CalcHiggsCorr(err)

	call setmudim(oldmudim)

	if( err .ne. 0 ) then
	  Warning("Could not evaluate corrections at mudim = 2 MT")
	  goto 999
	endif

	do i = 1, NHiggsErr
	  if( HiggsData(i) .ne. uncomputable .and.
     &        HiggsErr(i) .ne. uncomputable )
     &      HiggsErr(i) = max(HiggsErr(i),
     &        abs(HiggsData(i) - HiggsBase(i)))
	enddo

	do i = 1, NSfData
	  SfSave(i) = SfData(i)
	enddo

* b) exclude higher-order MB resummation effects

	if( debuglevel .ge. 2 )
     &    DHIGGS "Uncertainties: no higher-order MB resummation" ENDL

	if( botResum .eq. 0 ) then
	  bM = bTR
	  bM0 = bTR0
	  tl_valid = 0
	  call CalcHiggsCorr(err)
	  bM = 4
	  bM0 = 4
	  do i = 1, NHiggsErr
	    HiggsBase(i) = HiggsData(i)
	  enddo
	endif

	botResum = ibset(botResum, 8)
	tl_valid = 0
	call CalcHiggsCorr(err)
	botResum = ibclr(botResum, 8)

	if( err .ne. 0 ) then
	  Warning("Could not evaluate corrections for unresummed MB")
	  goto 998
	endif

	do i = 1, NHiggsErr
	  if( HiggsData(i) .ne. uncomputable .and.
     &        HiggsErr(i) .ne. uncomputable )
     &      HiggsErr(i) = HiggsErr(i) +
     &        abs(HiggsData(i) - HiggsBase(i))
	enddo

* c) use MTpole instead of MTrun in the two-loop corrections

	if( debuglevel .ge. 2 )
     &    DHIGGS "Uncertainties: MTpole in 2L corrections" ENDL

	runningMT = 1 - runningMT
	tM2 = 3 + (tT - 3)*runningMT
	if( looplevel .lt. 3 ) tM1 = tM2
	tl_valid = 0
	call CalcHiggsCorr(err)
	runningMT = 1 - runningMT
	tM1 = 3 + (tT - 3)*runningMT
	tM2 = tM1

	if( err .ne. 0 ) then
	  Warning("Could not evaluate corrections for MTrun = MTpole")
	  goto 998
	endif

	do i = 1, NHiggsErr
	  if( HiggsData(i) .ne. uncomputable .and.
     &        HiggsErr(i) .ne. uncomputable )
     &      HiggsErr(i) = HiggsErr(i) +
     &        abs(HiggsData(i) - HiggsBase(i))
	enddo

998	do i = 1, NSfData
	  SfData(i) = SfSave(i)
	enddo

999	do i = 1, NHiggsData
	  HiggsData(i) = HiggsSave(i)
	enddo
	tl_valid = valid
	higgs_valid = valid

	if( err .eq. 0 ) then
	  do i = 1, NHiggs
	    DeltaMHiggs_(i) = DeltaMHiggs(i)
	  enddo
	  DeltaSAeff_ = DeltaSAeff
	  do i = 1, NNeutral
	    do j = 1, NNeutral
	      DeltaUHiggs_(j,i) = DeltaXHiggs(j,i,1)
	      DeltaZHiggs_(j,i) = DeltaXHiggs(j,i,2)
	    enddo
	  enddo
	endif
	end
