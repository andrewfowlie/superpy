* TLhr.h
* declarations for Heidi Rzehak's two-loop corrections
* this file is part of FeynHiggs
* last modified 30 Apr 12 th

#include "FH.h"
#include "looptools.h"
#include "debug.h"
#include "TLhrvars.h"

#ifndef DTLHR
#define DTLHR if( debuglevel .ge. 2 ) DSELF

#define AtC Conjugate(At)
#define PhiAtC Conjugate(PhiAt)
#define MUE1C Conjugate(MUE1)
#define UStauC(i,j) Conjugate(UStau(i,j))
#define UStopC(i,j) Conjugate(UStop(i,j))
#define UCStopC(i,j) Conjugate(UCStop(i,j))
#define UUStopC(i,j) Conjugate(UUStop(i,j))
#endif

	RealType dMTfin
	parameter (dMTfin = 0)

	ComplexType Xt, At, PhiAt, MUEXt
	ComplexType UStop(2,2)
	ComplexType UCStop(3,4), UUStop(3,4)
	RealType UStop2(2,2)

	RealType MStop(2), MStop2(4)
	RealType MSbot(2), MSbot2(4)

	RealType MSq2Diff(2,2)
	RealType MTy, MTy2, MGlmT2, MGlpT2
	RealType MGlpTmSt2(2), MGlpTmSt4(2)
	RealType MGlpTmStxGlT4(2)
	RealType Q, MUE2
	RealType A0delStop(2), A0delGl, A0delT

	common /hrvar2/ Xt, At, PhiAt, MUEXt
	common /hrvar2/ UStop, UCStop, UUStop, UStop2
	common /hrvar2/ MStop, MStop2, MSbot, MSbot2, MSq2Diff
	common /hrvar2/ MTy, MTy2, MGlmT2, MGlpT2
	common /hrvar2/ MGlpTmSt2, MGlpTmSt4
	common /hrvar2/ MGlpTmStxGlT4, Q, MUE2
	common /hrvar2/ A0delStop, A0delGl, A0delT

	RealType T134
	external T134

