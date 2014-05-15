* gm2.h
* definitions for the g-2 calculation
* this file is part of FeynHiggs
* last modified 7 Jul 12 th

#include "FH.h"
#include "looptools.h"


	RealType Atop, Abot, Atau
	equivalence (Af(2,3), Atau)
	equivalence (Af(3,3), Atop)
	equivalence (Af(4,3), Abot)

	RealType MSneu(2), MSNeu2(4)
	equivalence (MSf(1,1,3), MSneu)
	equivalence (MSf2(1,1,3), MSneu2)

	RealType MStau(2), MStau2(4)
	ComplexType UStau11, UStau21, UStau12, UStau22
	equivalence (MSf(1,2,3), MStau)
	equivalence (MSf2(1,2,3), MStau2)
	equivalence (USf(1,1,2,3), UStau11)
	equivalence (USf(2,1,2,3), UStau21)
	equivalence (USf(1,2,2,3), UStau12)
	equivalence (USf(2,2,2,3), UStau22)

	RealType Mtop, Mtop2, MStop(2), MStop2(4)
	ComplexType UStop11, UStop21, UStop12, UStop22
	equivalence (MT, Mtop)
	equivalence (MT2, Mtop2)
	equivalence (MSf(1,3,3), MStop)
	equivalence (MSf2(1,3,3), MStop2)
	equivalence (USf(1,1,3,3), UStop11)
	equivalence (USf(2,1,3,3), UStop21)
	equivalence (USf(1,2,3,3), UStop12)
	equivalence (USf(2,2,3,3), UStop22)

	RealType Mbot, Mbot2, MSbot(2), MSbot2(4)
	ComplexType USbot11, USbot21, USbot12, USbot22
	equivalence (Mf(bBR,3), Mbot)
	equivalence (Mf2(bBR,3), Mbot2)
	equivalence (MSf(1,bBR,3), MSbot)
	equivalence (MSf2(1,bBR,3), MSbot2)
	equivalence (USf(1,1,bBR,3), USbot11)
	equivalence (USf(2,1,bBR,3), USbot21)
	equivalence (USf(1,2,bBR,3), USbot12)
	equivalence (USf(2,2,bBR,3), USbot22)

	ComplexType gm2_1L, gm2_2L
	RealType MSl2Diff(2,1), MSq2Diff(2,2)

	common /gm2_global/
     &    gm2_1L, gm2_2L,
     &    MSl2Diff, MSq2Diff

	ComplexType TF
	external TF

