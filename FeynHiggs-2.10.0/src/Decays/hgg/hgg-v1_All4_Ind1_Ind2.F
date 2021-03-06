#if 0
* hgg-v1_All4_Ind1_Ind2.F
* form factors for {S[10, {hno}]} -> {V[5, {Glu2}], V[5, {Glu3}]}
* generated by FormCalc 8.3 on 3-Dec-2013 16:50
* this file is part of FeynHiggs
#endif

#include "hgg-vars.h"
#include "inline.h"

	subroutine hgg_v1_All4_Ind1_Ind2
	implicit none

#include "hgg-vars.h"

#include "inline.h"

        Cloop(1) = 
     &   Cloop(1) - 
     &    1/(12.D0*Pi)*(AlfasMH*EL0*
     &        (Sub22(All4,Ind1,Ind2)*Sub3(All4) - 
     &          Sub25(All4)*Sub44(All4,Ind1,Ind2)))/(MW*SW)

#ifdef DETAILED_DEBUG
	DCOUP "v1_All4_Ind1_Ind2: Cloop(1) =", Cloop(1) ENDL
#endif

	end
