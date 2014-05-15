* util.h
* prototypes for the functions in util.a
* this file is part of FormCalc
* last modified 7 Apr 04 th


	integer Delta
	double precision ThreeMom, Li2
	double precision SInvariant, TInvariant
	double complex Pair, Eps
	double complex SxS, SeS
	integer VxS, VeS, BxS, BeS

	external Delta
	external ThreeMom, Li2
	external SInvariant, TInvariant
	external Pair, Eps
	external SxS, SeS
	external VxS, VeS, BxS, BeS

#ifdef LEGS
	double complex vec(2,2, 8, 0:LEGS)
	common /vectors/ vec

	double precision momspec(4, LEGS)
	common /momenta/ momspec
#endif

#ifndef k
#define k(i) (8*i+1)
#define s(i) (8*i+3)
#define e(i) (8*i+3+Hel(i))
#define ec(i) (8*i+3-Hel(i))
#define Spinor(i,s,om) (s*2*Hel(i)+16*i+om+5)
#define DottedSpinor(i,s,om) (s*2*Hel(i)+16*i+om+7)
#endif

