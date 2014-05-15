* Decay.h
* common definitions for all decays
* this file is part of FeynHiggs
* last modified 21 Sep 12 th

#include "FH.h"
#include "looptools.h"

#ifndef k
#define LEGS 3

#define k(i) (nvec*(i-1)+1)
#define s(i) (nvec*(i-1)+3)
#define e(i) (nvec*(i-1)+3+Hel(i))
#define ec(i) (nvec*(i-1)+3-Hel(i))
#define Spinor(i,s,d) (s*Hel(i)+nvec*(i-1)+d+5)

#define TEST(i,b) if( btest(i,b) ) then
#define ENDTEST(i,b) endif

#define BIT_HEL(i) (3*(LEGS-i)+Hel(i)+1)
#define LOOP_HEL(h) do h = -1, 1
#define ENDLOOP_HEL(h) enddo
#endif

	integer nvec
	parameter (nvec = 8)

	ComplexType vec(2,2,nvec*LEGS)
	common /vectors/ vec

	RealType mass(2,LEGS)
	common /masses/ mass

	RealType m1, m2, m3, m12, m22, m32
	equivalence (mass(1,1), m1), (mass(2,1), m12)
	equivalence (mass(1,2), m2), (mass(2,2), m22)
	equivalence (mass(1,3), m3), (mass(2,3), m32)

	ComplexType HffDb(0:1,3,2:4,3)
	RealType AlfasMH, Divergence
	integer hno, hno1, hno2, gno1, gno2, sub1L
	common /decay/ HffDb
	common /decay/ AlfasMH, Divergence
	common /decay/ hno, hno1, hno2, gno1, gno2, sub1L
