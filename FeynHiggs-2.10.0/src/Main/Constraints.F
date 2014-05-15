* Constraints.F
* evaluation of electroweak precision observables
* used as further constraints for the MSSM parameter space
* this file is part of FeynHiggs
* last modified 26 Sep 12 th

#include "externals.h"
#include "types.h"
#include "debug.h"


	subroutine FHConstraints(error, gm2,
     &    deltarho, MWMSSM, MWSM, SW2MSSM, SW2SM,
     &    edmeTh, edmn, edmHg,
     &    ccb)
	implicit none
	integer error
	RealType gm2
	RealType deltarho, MWMSSM, MWSM, SW2MSSM, SW2SM
	RealType edmeTh, edmn, edmHg
	integer ccb

#include "gm2.h"
#define __SUBROUTINE__ "FHConstraints"

	RealType LambdaNP, gm2_2L_QEDfactor
	RealType cu, cd

	call CheckCouplings(error, 1)

	if( error .gt. 0 ) return

* a) g-2, see hep-ph/0312264

* two-loop corrections from QED logs, hep-ph/9803384
	LambdaNP = (MSneu(1) + MStau(1) + MStau(2) + 
     &    MCha(1) + MCha(2) +
     &    MNeu(1) + MNeu(2) + MNeu(3) + MNeu(4))/9D0
	gm2_2L_QEDfactor = -4*AlfaGF/pi*log(LambdaNP/MM)

* two-loop corrections from hep-ph/0312264
	MSl2Diff(1,1) = (MStau(1) - MSneu(1))*(MStau(1) + MSneu(1))
	MSl2Diff(2,1) = (MStau(2) - MSneu(1))*(MStau(2) + MSneu(1))
	MSq2Diff(1,1) = (MSbot(1) - MStop(1))*(MSbot(1) + MStop(1))
	MSq2Diff(2,1) = (MSbot(2) - MStop(1))*(MSbot(2) + MStop(1))
	MSq2Diff(1,2) = (MSbot(1) - MStop(2))*(MSbot(1) + MStop(2))
	MSq2Diff(2,2) = (MSbot(2) - MStop(2))*(MSbot(2) + MStop(2))

	if( min(abs(MNeu(1) - MSf(1,2,2)),
     &          abs(MNeu(2) - MSf(1,2,2)),
     &          abs(MNeu(3) - MSf(1,2,2)),
     &          abs(MNeu(4) - MSf(1,2,2)),
     &          abs(MCha(1) - MSf(1,2,2)),
     &          abs(MCha(2) - MSf(1,2,2)),
     &          abs(MNeu(1) - MSf(2,2,2)),
     &          abs(MNeu(2) - MSf(2,2,2)),
     &          abs(MNeu(3) - MSf(2,2,2)),
     &          abs(MNeu(4) - MSf(2,2,2)),
     &          abs(MCha(1) - MSf(2,2,2)),
     &          abs(MCha(2) - MSf(2,2,2))) .lt. .1D0 )
     &    Warning("g-2 possibly unreliable due to numerical problems")

	call gm2_1L_1

	call gm2_2L_1
	call gm2_2L_2
	if( debuglevel .ge. 2 )
     &    DCONST "gm2 2-loop short =", gm2_2L ENDL

#if GM2FULL
	call gm2_2L_full1
	call gm2_2L_full2
	call gm2_2L_full3
	call gm2_2L_full4
	call gm2_2L_full5
	call gm2_2L_full6
	call gm2_2L_full7
	call gm2_2L_full8
	call gm2_2L_full9
	call gm2_2L_full10
	call gm2_2L_full11
	if( debuglevel .ge. 2 )
     &    DCONST "gm2 2-loop full  =", gm2_2L ENDL
#endif

	gm2 = (1 + gm2_2L_QEDfactor)*Re(gm2_1L) +
     &    5/(24*pi**2)*AlfaGF**2*MM2/(MW2*SW2)*Re(gm2_2L)

* b) EW precision observables
	call EWPO(deltarho, MWMSSM, MWSM, SW2MSSM, SW2SM)

* c) EDMs
	call CalcEDMs(edmeTh, edmn, edmHg)

* d) test for colour-breaking minimum
	cd = Sq(MUE) - .5D0*MZ2
	cu = (MA02 + MZ2)*CB2 + cd
	cd = (MA02 + MZ2)*SB2 + cd
	ccb = 0
	if( Sq(Af(2,3)) .gt. 3*(Re(MSS2(3,3,1)) +
     &                          Re(MSS2(3,3,2)) + cd) .or.
     &      Sq(Af(3,3)) .gt. 3*(Re(MSS2(3,3,3)) +
     &                          Re(MSS2(3,3,4)) + cu) .or.
     &      Sq(Af(4,3)) .gt. 3*(Re(MSS2(3,3,3)) +
     &                          Re(MSS2(3,3,5)) + cd) )
     &    ccb = 1
	end

