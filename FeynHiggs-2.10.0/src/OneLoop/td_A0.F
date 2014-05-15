* td_A0.F
* this file is part of FeynHiggs
* generated 05 Nov 2013 19:22

#include "externals.h"
#include "types.h"
#include "debug.h"


	ComplexType function td_A0(bX)
	implicit none
	integer bX


#include "FH.h"

	ComplexType td_A0_mfv, td_A0_nmfv
	external td_A0_mfv, td_A0_nmfv

	if( ibits(fv, 3, 2) .eq. 0 ) then
	  td_A0 = td_A0_mfv(bX)
	else
	  td_A0 = td_A0_nmfv(bX)
	endif

#ifdef DEBUG
	DHIGGS 'td_A0 =', td_A0 ENDL
#endif
	end


************************************************************************


	ComplexType function td_A0_mfv(bX)
	implicit none
	integer bX

#include "FH.h"
#include "looptools.h"

	integer Cha2, Gen2, Neu2, Sfe2


	td_A0_mfv = 0

	LOOP(Gen2, gM,3,1)
	LOOP(Sfe2, 1,2,1)

        td_A0_mfv = td_A0_mfv + 
     &    (3*cI)/(32.D0*Pi**2)*
     &     (EL1L*A0(MSf2(Sfe2,tM1,Gen2))*
     &        ((Kf(Gen2,Gen2,tM1) + MUEC*TB*Mf(tM1,Gen2))*
     &           USf(Sfe2,2,tM1,Gen2)*USfC(Sfe2,1,tM1,Gen2) - 
     &          (KfC(Gen2,Gen2,tM1) + MUE*TB*Mf(tM1,Gen2))*
     &           USf(Sfe2,1,tM1,Gen2)*USfC(Sfe2,2,tM1,Gen2)))/
     &      (MW*SW*TB)

	ENDLOOP(Sfe2)
	ENDLOOP(Gen2)

#ifdef DETAILED_DEBUG
	DHIGGS "td_A0_mfv =", td_A0_mfv ENDL
#endif

	if( mssmpart .le. 1 ) return

	LOOP(Gen2, gM,3,1)
	LOOP(Sfe2, 1,2,1)

        td_A0_mfv = td_A0_mfv + 
     &    (3*cI)/(32.D0*Pi**2)*
     &     (EL1L*A0(MSf2(Sfe2,bX,Gen2))*
     &        ((TB*Kf(Gen2,Gen2,bX) + MUEC*Mf(bX,Gen2))*
     &           USf(Sfe2,2,bX,Gen2)*USfC(Sfe2,1,bX,Gen2) - 
     &          (TB*KfC(Gen2,Gen2,bX) + MUE*Mf(bX,Gen2))*
     &           USf(Sfe2,1,bX,Gen2)*USfC(Sfe2,2,bX,Gen2)))/(MW*SW)

	ENDLOOP(Sfe2)
	ENDLOOP(Gen2)

#ifdef DETAILED_DEBUG
	DHIGGS "td_A0_mfv =", td_A0_mfv ENDL
#endif

	if( mssmpart .le. 2 ) return

	LOOP(Sfe2, 1,2,1)
	LOOP(Gen2, 1,3,1)

        td_A0_mfv = td_A0_mfv + 
     &    cI/(32.D0*Pi**2)*(EL1L*A0(MSf2(Sfe2,2,Gen2))*
     &        ((TB*Kf(Gen2,Gen2,2) + MUEC*Mf(2,Gen2))*
     &           USf(Sfe2,2,2,Gen2)*USfC(Sfe2,1,2,Gen2) - 
     &          (TB*KfC(Gen2,Gen2,2) + MUE*Mf(2,Gen2))*
     &           USf(Sfe2,1,2,Gen2)*USfC(Sfe2,2,2,Gen2)))/(MW*SW)

	ENDLOOP(Gen2)
	ENDLOOP(Sfe2)

#ifdef DETAILED_DEBUG
	DHIGGS "td_A0_mfv =", td_A0_mfv ENDL
#endif

	if( mssmpart .le. 3 ) return

	LOOP(Cha2, 1,2,1)

        td_A0_mfv = td_A0_mfv - 
     &    cI/(8.D0*Pi**2*sqrt2)*
     &     (EL1L*A0(MCha2(Cha2))*MCha(Cha2)*
     &        (SB*(UCha(Cha2,2)*VCha(Cha2,1) - 
     &             UChaC(Cha2,2)*VChaC(Cha2,1)) + 
     &          CB*(UCha(Cha2,1)*VCha(Cha2,2) - 
     &             UChaC(Cha2,1)*VChaC(Cha2,2))))/SW

	ENDLOOP(Cha2)

	LOOP(Neu2, 1,4,1)

        td_A0_mfv = td_A0_mfv + 
     &    cI/(16.D0*Pi**2)*(EL1L*A0(MNeu2(Neu2))*MNeu(Neu2)*
     &        ((SW*ZNeu(Neu2,1) - CW*ZNeu(Neu2,2))*
     &           (SB*ZNeu(Neu2,3) - CB*ZNeu(Neu2,4)) - 
     &          (SW*ZNeuC(Neu2,1) - CW*ZNeuC(Neu2,2))*
     &           (SB*ZNeuC(Neu2,3) - CB*ZNeuC(Neu2,4))))/(CW*SW)

	ENDLOOP(Neu2)

#ifdef DETAILED_DEBUG
	DHIGGS "td_A0_mfv =", td_A0_mfv ENDL
#endif

	end


************************************************************************


	ComplexType function td_A0_nmfv(bX)
	implicit none
	integer bX

#include "FH.h"
#include "looptools.h"

	integer All2, Cha2, Ind1, Ind2, Neu2

	ComplexType tmp1, tmp2, tmp3

	td_A0_nmfv = 0

	LOOP(All2, 1,6,gM)

	tmp1 = A0(MASf2(All2,tM1))

	LOOP(Ind2, 1,3,1)
	LOOP(Ind1, 1,3,1)

        td_A0_nmfv = td_A0_nmfv + 
     &    (3*cI)/(32.D0*Pi**2)*
     &     (EL1L*tmp1*((Kf(Ind1,Ind2,tM1) + 
     &             MUEC*TB*Delta(Ind1,Ind2)*Mf(tM1,Ind1))*
     &           UASf(All2,3 + Ind2,tM1)*UASfC(All2,Ind1,tM1) - 
     &          (KfC(Ind2,Ind1,tM1) + 
     &             MUE*TB*Delta(Ind1,Ind2)*Mf(tM1,Ind1))*
     &           UASf(All2,Ind2,tM1)*UASfC(All2,3 + Ind1,tM1)))/
     &      (MW*SW*TB)

	ENDLOOP(Ind1)
	ENDLOOP(Ind2)

	ENDLOOP(All2)

#ifdef DETAILED_DEBUG
	DHIGGS "td_A0_nmfv =", td_A0_nmfv ENDL
#endif

	if( mssmpart .le. 1 ) return

	LOOP(All2, 1,6,gM)

	tmp2 = A0(MASf2(All2,bX))

	LOOP(Ind2, 1,3,1)
	LOOP(Ind1, 1,3,1)

        td_A0_nmfv = td_A0_nmfv + 
     &    (3*cI)/(32.D0*Pi**2)*
     &     (EL1L*tmp2*((TB*Kf(Ind1,Ind2,bX) + 
     &             MUEC*Delta(Ind1,Ind2)*Mf(bX,Ind1))*
     &           UASf(All2,3 + Ind2,bX)*UASfC(All2,Ind1,bX) - 
     &          (TB*KfC(Ind2,Ind1,bX) + 
     &             MUE*Delta(Ind1,Ind2)*Mf(bX,Ind1))*
     &           UASf(All2,Ind2,bX)*UASfC(All2,3 + Ind1,bX)))/
     &      (MW*SW)

	ENDLOOP(Ind1)
	ENDLOOP(Ind2)

	ENDLOOP(All2)

#ifdef DETAILED_DEBUG
	DHIGGS "td_A0_nmfv =", td_A0_nmfv ENDL
#endif

	if( mssmpart .le. 2 ) return

	LOOP(All2, 1,6,1)

	tmp3 = A0(MASf2(All2,2))

	LOOP(Ind2, 1,3,1)
	LOOP(Ind1, 1,3,1)

        td_A0_nmfv = td_A0_nmfv + 
     &    cI/(32.D0*Pi**2)*(EL1L*tmp3*
     &        ((TB*Kf(Ind1,Ind2,2) + 
     &             MUEC*Delta(Ind1,Ind2)*Mf(2,Ind1))*
     &           UASf(All2,3 + Ind2,2)*UASfC(All2,Ind1,2) - 
     &          (TB*KfC(Ind2,Ind1,2) + 
     &             MUE*Delta(Ind1,Ind2)*Mf(2,Ind1))*
     &           UASf(All2,Ind2,2)*UASfC(All2,3 + Ind1,2)))/(MW*SW)

	ENDLOOP(Ind1)
	ENDLOOP(Ind2)

	ENDLOOP(All2)

#ifdef DETAILED_DEBUG
	DHIGGS "td_A0_nmfv =", td_A0_nmfv ENDL
#endif

	if( mssmpart .le. 3 ) return

	LOOP(Cha2, 1,2,1)

        td_A0_nmfv = td_A0_nmfv - 
     &    cI/(8.D0*Pi**2*sqrt2)*
     &     (EL1L*A0(MCha2(Cha2))*MCha(Cha2)*
     &        (SB*(UCha(Cha2,2)*VCha(Cha2,1) - 
     &             UChaC(Cha2,2)*VChaC(Cha2,1)) + 
     &          CB*(UCha(Cha2,1)*VCha(Cha2,2) - 
     &             UChaC(Cha2,1)*VChaC(Cha2,2))))/SW

	ENDLOOP(Cha2)

	LOOP(Neu2, 1,4,1)

        td_A0_nmfv = td_A0_nmfv + 
     &    cI/(16.D0*Pi**2)*(EL1L*A0(MNeu2(Neu2))*MNeu(Neu2)*
     &        ((SW*ZNeu(Neu2,1) - CW*ZNeu(Neu2,2))*
     &           (SB*ZNeu(Neu2,3) - CB*ZNeu(Neu2,4)) - 
     &          (SW*ZNeuC(Neu2,1) - CW*ZNeuC(Neu2,2))*
     &           (SB*ZNeuC(Neu2,3) - CB*ZNeuC(Neu2,4))))/(CW*SW)

	ENDLOOP(Neu2)

#ifdef DETAILED_DEBUG
	DHIGGS "td_A0_nmfv =", td_A0_nmfv ENDL
#endif

	end

