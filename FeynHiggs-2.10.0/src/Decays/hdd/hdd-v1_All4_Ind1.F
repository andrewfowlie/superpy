#if 0
* hdd-v1_All4_Ind1.F
* form factors for {S[10, {hno}]} -> {F[4, {gno1, Col2}], -F[4, {gno2, Col3}]}
* generated by FormCalc 8.3 on 3-Dec-2013 17:01
* this file is part of FeynHiggs
#endif

#include "hdd-vars.h"
#include "inline.h"

	subroutine hdd_v1_All4_Ind1
	implicit none

#include "hdd-vars.h"

#include "inline.h"

        Cloop(1) = 
     &   Cloop(1) + 
     &    cI/(32.D0*Pi)*(AlfaMZ*ELMZ*(F1 - F2)*S2B*Delta(gno1,gno2)*
     &        Mf(bH,gno1)*(SB2*
     &           (pave7(All4)*Sub113(All4,Ind1) + 
     &             pave8(All4)*Sub116(All4,Ind1)) - 
     &          CB2*pave9(All4)*Sub119(All4,Ind1))*
     &        XHiggs(hno,3,uzext))/
     &      (CB2*CW2*MW*MW2*(m12 - MZ2)*SB2*SW*SW2)

#ifdef DETAILED_DEBUG
	DCOUP "v1_All4_Ind1: Cloop(1) =", Cloop(1) ENDL
#endif

	end
