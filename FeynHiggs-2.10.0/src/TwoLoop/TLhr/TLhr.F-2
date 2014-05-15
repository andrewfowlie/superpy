* TLhr.F
* master file for Heidi Rzehak's two-loop corrections
* this file is part of FeynHiggs
* last modified 13 Nov 13 th

#include "externals.h"
#include "types.h"


	subroutine TLhrSE(error, mudim)
	implicit none
	integer error
	RealType mudim

#include "TLhr.h"
#define __SUBROUTINE__ "TLhrSE"

	ComplexType seinput
	RealType c

	call TLMSf_gaugeless(error, MUE, Xf(3,3))

	Q = mudim
	c = abs(At)
	PhiAt = c
	if( c .gt. 1D-14 ) PhiAt = At/c
	MUE2 = Sq(MUE)

	UUStop(1,1) = UStop(1,1)**2
	UUStop(2,1) = UStop(2,1)**2
	UUStop(1,2) = UStop(1,2)**2
	UUStop(2,2) = UStop(2,2)**2

	UUStop(3,1) = UStop(1,1)*UStop(2,1)
	UUStop(3,2) = UStop(1,2)*UStop(2,2)
	UUStop(1,3) = UStop(1,1)*UStop(1,2)
	UUStop(2,3) = UStop(2,1)*UStop(2,2)
	UUStop(3,3) = UStop(1,1)*UStop(2,2)
	UUStop(3,4) = UStop(1,2)*UStop(2,1)

	c = 2*MGl*MTy
	MGlmT2 = MGl2 - MTy2
	MGlpT2 = MGl2 + MTy2
	MGlpTmSt2(1) = MGlpT2 - MStop2(1)
	MGlpTmSt4(1) = MGlpTmSt2(1)**2
	MGlpTmStxGlT4(1) = (MGlpTmSt2(1) - c)*(MGlpTmSt2(1) + c)
	MGlpTmSt2(2) = MGlpT2 - MStop2(2)
	MGlpTmSt4(2) = MGlpTmSt2(2)**2
	MGlpTmStxGlT4(2) = (MGlpTmSt2(2) - c)*(MGlpTmSt2(2) + c)

	A0delStop(1) = A0delq(MStop2(1), Q)
	A0delStop(2) = A0delq(MStop2(2), Q)
	A0delGl = A0delq(MGl2, Q)
	A0delT = A0delq(MTy2, Q)

	if( debuglevel .ge. 2 ) then
	  DSELF "==========================="	ENDL
	  DSELF "two-loop self-energies (hr)"	ENDL
	  DSELF "---------------------------"	ENDL
	endif

	call TLhr_abbr0
	call TLhr_rcs

	call TLhr_abbr1_1
	call TLhr_abbr1_2
	call TLhr_abbr1_3
	call TLhr_abbr1_4
	call TLhr_abbr1_5
	call TLhr_abbr1_6
	call TLhr_abbr1_7
	call TLhr_abbr1_8
	call TLhr_abbr1_9
	call TLhr_abbr1_10
	call TLhr_abbr1_11
	call TLhr_abbr1_12
	call TLhr_abbr1_13
	call TLhr_abbr1_14
	call TLhr_abbr1_15
	call TLhr_abbr1_16
	call TLhr_abbr1_17
	call TLhr_abbr1_18
	call TLhr_self1

	c = EL2L/(2*SW*MW)

	if( inputmass .eq. 3 ) then
	  seinput = selfA0A0
	else
	  seinput = selfHmHp
	endif

	se2R(h0h0) = selfh0h0 - CBA2*seinput +
     &    c*SBA*(CBA*SBA*tadHH - (1 + CBA2)*tadh0)

	se2R(HHHH) = selfHHHH - SBA2*seinput +
     &    c*CBA*(CBA*SBA*tadh0 - (1 + SBA2)*tadHH)

	se2R(A0A0) = selfA0A0 - seinput

	se2R(HmHp) = selfHmHp - seinput

	if( debuglevel .ge. 2 ) then
	  DSELF "se2R(h0h0) =", se2R(h0h0)	ENDL
	  DSELF "se2R(HHHH) =", se2R(HHHH)	ENDL
	  DSELF "se2R(A0A0) =", se2R(A0A0)	ENDL
	  DSELF "se2R(HmHp) =", se2R(HmHp)	ENDL
	endif

	if( higgsmix .lt. 2 ) return

	call TLhr_abbr2_1
	call TLhr_abbr2_2
	call TLhr_abbr2_3
	call TLhr_self2

	se2R(h0HH) = selfh0HH + CBA*SBA*seinput -
     &    c*(CBA**3*tadh0 + SBA**3*tadHH)

	if( debuglevel .ge. 2 ) then
	  DSELF "se2R(h0HH) =", se2R(h0HH)	ENDL
	endif

	if( higgsmix .lt. 3 ) return

	call TLhr_abbr3_1
	call TLhr_abbr3_2
	call TLhr_abbr3_3
	call TLhr_abbr3_4
	call TLhr_self3

	se2R(h0A0) = selfh0A0 - c*SBA*tadA0

	se2R(HHA0) = selfHHA0 - c*CBA*tadA0

	if( debuglevel .ge. 2 ) then
	  DSELF "se2R(h0A0) =", se2R(h0A0)	ENDL
	  DSELF "se2R(HHA0) =", se2R(HHA0)	ENDL
	endif
	end


************************************************************************

	subroutine TLMSf_gaugeless(error, MUE_, Xt_)
	implicit none
	integer error
	ComplexType MUE_, Xt_

#include "TLhr.h"

	RealType DSStop2(2)
	integer fail

	Xt = Xt_
	At = Xt_ + Conjugate(MUE_)/TB
	MUEXt = MUE_*Xt

	MTy = Mf(tM2,3)
	MTy2 = MTy**2

* use squark masses w/o D-terms & neglect MB
	DSStop2(1) = DSS2(1,3,3) - (DSf(1,3) - DSf(2,3))
	DSStop2(2) = DSS2(2,3,3) - DSf(2,3)
	if( DSStop2(1)**2 + Sq(Xf(3,3)) .lt. .1D0 )
     &    DSStop2(1) = -.2D0*sqrt(DSStop2(2))
	fail = 0
	call EigenSf(fail, 1, MStop, MStop2, UStop,
     &    DSStop2, Mf(tM1,3), Xf(3,3), 0)
	if( fail .ne. 0 ) then
	  Error(error, "negative stop mass squares")
	  return
	endif

	MSbot2(1) = DSS2(1,4,3) + DSS2(2,4,3) - DSf(1,4)
	MSbot2(2) = DSS2(2,4,3) - DSf(2,4)
	MSbot2(3) = MSbot2(2) - MSbot2(1)
	MSbot2(4) = MSbot2(2) + MSbot2(1)
	MSbot(1) = sqrt(MSbot2(1))
	MSbot(2) = sqrt(MSbot2(2))

	if( debuglevel .ge. 4 ) then
	  DHIGGS "gaugeless approximation:"	ENDL
	  DHIGGS "MStop =", MStop(1), MStop(2)	ENDL
	  DHIGGS "MSbot =", MSbot(1), MSbot(2)	ENDL
	endif

c	USbot(1,1) = 1
c	USbot(2,1) = 0
c	USbot(1,2) = 0
c	USbot(2,2) = 1

	MSq2Diff(1,1) = (MSbot(1) - MStop(1))*(MSbot(1) + MStop(1))
	MSq2Diff(2,1) = (MSbot(2) - MStop(1))*(MSbot(2) + MStop(1))
	MSq2Diff(1,2) = (MSbot(1) - MStop(2))*(MSbot(1) + MStop(2))
	MSq2Diff(2,2) = (MSbot(2) - MStop(2))*(MSbot(2) + MStop(2))

	UStop2(1,1) = Sq(UStop(1,1))
	UStop2(2,1) = Sq(UStop(2,1))
	UStop2(1,2) = Sq(UStop(1,2))
	UStop2(2,2) = Sq(UStop(2,2))

	UCStop(3,1) = UStop(1,1)*UStopC(2,1)
	UCStop(3,2) = UStop(1,2)*UStopC(2,2)
	UCStop(1,3) = UStop(1,1)*UStopC(1,2)
	UCStop(2,3) = UStop(2,1)*UStopC(2,2)
	UCStop(3,3) = UStop(1,1)*UStopC(2,2)
	UCStop(3,4) = UStop(1,2)*UStopC(2,1)
	end


