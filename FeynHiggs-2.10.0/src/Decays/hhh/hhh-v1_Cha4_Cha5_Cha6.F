#if 0
* hhh-v1_Cha4_Cha5_Cha6.F
* form factors for {S[10, {hno}]} -> {S[10, {hno1}], S[10, {hno2}]}
* generated by FormCalc 8.3 on 3-Dec-2013 16:52
* this file is part of FeynHiggs
#endif

#include "hhh-vars.h"
#include "inline.h"

	subroutine hhh_v1_Cha4_Cha5_Cha6
	implicit none

#include "hhh-vars.h"

#include "inline.h"

        Cloop(1) = 
     &   Cloop(1) + 
     &    cI/(8.D0*Pi*sqrt2)*
     &     (Alfa0*EL0*(cI*(pave11(cc0,Cha4,Cha5,Cha6)*
     &             Sub958(Cha4,Cha5,Cha6)) + 
     &          2*pave10(Cha5,Cha6)*Sub949(Cha4,Cha5,Cha6) - 
     &          Sub954(Cha4,Cha5,Cha6)))/(SW*SW2)

#ifdef DETAILED_DEBUG
	DCOUP "v1_Cha4_Cha5_Cha6: Cloop(1) =", Cloop(1) ENDL
#endif

	end
