#if 0
* hgaZ-v1.F
* form factors for {S[10, {hno}]} -> {V[1], V[2]}
* generated by FormCalc 8.3 on 3-Dec-2013 16:50
* this file is part of FeynHiggs
#endif

#include "hgaZ-vars.h"
#include "inline.h"

	subroutine hgaZ_v1
	implicit none

#include "hgaZ-vars.h"

#include "inline.h"

        Cloop(1) = 
     &   Cloop(1) + 
     &    1/(8.D0*Pi)*(EL0*(Alfa0*Sub184 + Pi*(CW*MW2*Sub193*SW2)))/
     &      (CW*MW*SW2)

#ifdef DETAILED_DEBUG
	DCOUP "v1: Cloop(1) =", Cloop(1) ENDL
#endif

	end
