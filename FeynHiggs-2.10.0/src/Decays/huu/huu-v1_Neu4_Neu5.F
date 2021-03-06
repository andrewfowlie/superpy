#if 0
* huu-v1_Neu4_Neu5.F
* form factors for {S[10, {hno}]} -> {F[3, {gno1, Col2}], -F[3, {gno2, Col3}]}
* generated by FormCalc 8.3 on 3-Dec-2013 16:57
* this file is part of FeynHiggs
#endif

#include "huu-vars.h"
#include "inline.h"

	subroutine huu_v1_Neu4_Neu5
	implicit none

#include "huu-vars.h"

#include "inline.h"

        Cloop(1) = 
     &   Cloop(1) + 
     &    cI/(32.D0*Pi)*(AlfaMZ*ELMZ*(F1 - F2)*Delta(gno1,gno2)*
     &        Mf(tH,gno1)*(CW*pave14(Neu5)*Sub311(Neu4,Neu5) + 
     &          pave15(Neu4,Neu5)*Sub316(Neu4,Neu5) + 
     &          pave16(Neu4,Neu5)*Sub319(Neu4,Neu5)))/
     &      (CW*CW2*MW*(m12 - MZ2)*SW*SW2)

#ifdef DETAILED_DEBUG
	DCOUP "v1_Neu4_Neu5: Cloop(1) =", Cloop(1) ENDL
#endif

	end
