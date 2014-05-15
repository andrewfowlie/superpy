* td_h0.F
* this file is part of FeynHiggs
* generated 05 Nov 2013 19:22

#include "externals.h"
#include "types.h"
#include "debug.h"


	ComplexType function td_h0(bX)
	implicit none
	integer bX


#include "FH.h"

	ComplexType td_h0_mfv, td_h0_nmfv
	external td_h0_mfv, td_h0_nmfv

	if( ibits(fv, 3, 2) .eq. 0 ) then
	  td_h0 = td_h0_mfv(bX)
	else
	  td_h0 = td_h0_nmfv(bX)
	endif

#ifdef DEBUG
	DHIGGS 'td_h0 =', td_h0 ENDL
#endif
	end


************************************************************************


	ComplexType function td_h0_mfv(bX)
	implicit none
	integer bX

#include "FH.h"
#include "looptools.h"

	integer Cha2, Gen2, Neu2, Sfe2


	td_h0_mfv = 0

	LOOP(Gen2, gM,3,1)

        td_h0_mfv = td_h0_mfv - 
     &    3/(8.D0*Pi**2)*(CA*EL1L*A0(Mf2(tM1,Gen2))*Mf2(tM1,Gen2))/
     &      (MW*SB*SW)

	ENDLOOP(Gen2)

	LOOP(Gen2, gM,3,1)
	LOOP(Sfe2, 1,2,1)

        td_h0_mfv = td_h0_mfv - 
     &    1/(32.D0*Pi**2)*(EL1L*A0(MSf2(Sfe2,tM1,Gen2))*
     &        (USf(Sfe2,1,tM1,Gen2)*
     &           ((MW*MZ*SAB*SB*(3 - 4*SW2) - 
     &                6*CA*CW*Mf2(tM1,Gen2))*USfC(Sfe2,1,tM1,Gen2)-
     &               3*CW*
     &              (CA*KfC(Gen2,Gen2,tM1) + MUE*SA*Mf(tM1,Gen2))*
     &              USfC(Sfe2,2,tM1,Gen2)) - 
     &          USf(Sfe2,2,tM1,Gen2)*
     &           (3*CW*(CA*Kf(Gen2,Gen2,tM1) + 
     &                MUEC*SA*Mf(tM1,Gen2))*USfC(Sfe2,1,tM1,Gen2)-
     &               2*(2*MW*MZ*SAB*SB*SW2 - 
     &                3*CA*CW*Mf2(tM1,Gen2))*USfC(Sfe2,2,tM1,Gen2))
     &          ))/(CW*MW*SB*SW)

	ENDLOOP(Sfe2)
	ENDLOOP(Gen2)

#ifdef DETAILED_DEBUG
	DHIGGS "td_h0_mfv =", td_h0_mfv ENDL
#endif

	if( mssmpart .le. 1 ) return

	LOOP(Gen2, gM,3,1)

        td_h0_mfv = td_h0_mfv + 
     &    3/(8.D0*Pi**2)*(EL1L*SA*A0(Mf2(bX,Gen2))*Mf2(bX,Gen2))/
     &      (CB*MW*SW)

	ENDLOOP(Gen2)

	LOOP(Gen2, gM,3,1)
	LOOP(Sfe2, 1,2,1)

        td_h0_mfv = td_h0_mfv + 
     &    1/(32.D0*Pi**2)*(EL1L*A0(MSf2(Sfe2,bX,Gen2))*
     &        (USf(Sfe2,1,bX,Gen2)*
     &           ((CB*MW*MZ*SAB*(3 - 2*SW2) - 
     &                6*CW*SA*Mf2(bX,Gen2))*USfC(Sfe2,1,bX,Gen2) - 
     &             3*CW*(SA*KfC(Gen2,Gen2,bX) + 
     &                CA*MUE*Mf(bX,Gen2))*USfC(Sfe2,2,bX,Gen2)) - 
     &          USf(Sfe2,2,bX,Gen2)*
     &           (3*CW*(SA*Kf(Gen2,Gen2,bX) + CA*MUEC*Mf(bX,Gen2))*
     &              USfC(Sfe2,1,bX,Gen2) - 
     &             2*(CB*MW*MZ*SAB*SW2 - 3*CW*SA*Mf2(bX,Gen2))*
     &              USfC(Sfe2,2,bX,Gen2))))/(CB*CW*MW*SW)

	ENDLOOP(Sfe2)
	ENDLOOP(Gen2)

#ifdef DETAILED_DEBUG
	DHIGGS "td_h0_mfv =", td_h0_mfv ENDL
#endif

	if( mssmpart .le. 2 ) return

	LOOP(Gen2, 1,3,1)

        td_h0_mfv = td_h0_mfv - 
     &    1/(32.D0*Pi**2)*(EL1L*
     &        (CB*MW*MZ*SAB*A0(MSf2(1,1,Gen2)) - 
     &          4*CW*SA*A0(Mf2(2,Gen2))*Mf2(2,Gen2)))/(CB*CW*MW*SW)

	ENDLOOP(Gen2)

	LOOP(Sfe2, 1,2,1)
	LOOP(Gen2, 1,3,1)

        td_h0_mfv = td_h0_mfv + 
     &    1/(32.D0*Pi**2)*(EL1L*A0(MSf2(Sfe2,2,Gen2))*
     &        (USf(Sfe2,1,2,Gen2)*
     &           ((CB*MW*MZ*SAB*(1 - 2*SW2) - 2*CW*SA*Mf2(2,Gen2))*
     &              USfC(Sfe2,1,2,Gen2) - 
     &             CW*(SA*KfC(Gen2,Gen2,2) + CA*MUE*Mf(2,Gen2))*
     &              USfC(Sfe2,2,2,Gen2)) - 
     &          USf(Sfe2,2,2,Gen2)*
     &           (CW*(SA*Kf(Gen2,Gen2,2) + CA*MUEC*Mf(2,Gen2))*
     &              USfC(Sfe2,1,2,Gen2) - 
     &             2*(CB*MW*MZ*SAB*SW2 - CW*SA*Mf2(2,Gen2))*
     &              USfC(Sfe2,2,2,Gen2))))/(CB*CW*MW*SW)

	ENDLOOP(Gen2)
	ENDLOOP(Sfe2)

#ifdef DETAILED_DEBUG
	DHIGGS "td_h0_mfv =", td_h0_mfv ENDL
#endif

	if( mssmpart .le. 3 ) return

        td_h0_mfv = td_h0_mfv + 
     &    1/(64.D0*Pi**2)*(EL1L*MW*
     &        (SAB*(C2B*A0(MA02) + 3*C2A*A0(Mh02)) - 
     &          (2*CAB*S2A + C2A*SAB)*A0(MHH2) + 
     &          2*(C2B*SAB + 2*CW2*SBA)*A0(MHp2) - 
     &          2*(C2B*SAB - 6*CW2*SBA)*A0(MW2) - 
     &          (C2B*SAB - 6*SBA)*A0(MZ2)))/(CW2*SW)

	LOOP(Cha2, 1,2,1)

        td_h0_mfv = td_h0_mfv + 
     &    1/(8.D0*Pi**2*sqrt2)*
     &     (EL1L*A0(MCha2(Cha2))*MCha(Cha2)*
     &        (SA*(UCha(Cha2,2)*VCha(Cha2,1) + 
     &             UChaC(Cha2,2)*VChaC(Cha2,1)) - 
     &          CA*(UCha(Cha2,1)*VCha(Cha2,2) + 
     &             UChaC(Cha2,1)*VChaC(Cha2,2))))/SW

	ENDLOOP(Cha2)

	LOOP(Neu2, 1,4,1)

        td_h0_mfv = td_h0_mfv - 
     &    1/(16.D0*Pi**2)*(EL1L*A0(MNeu2(Neu2))*MNeu(Neu2)*
     &        ((SW*ZNeu(Neu2,1) - CW*ZNeu(Neu2,2))*
     &           (SA*ZNeu(Neu2,3) + CA*ZNeu(Neu2,4)) + 
     &          (SW*ZNeuC(Neu2,1) - CW*ZNeuC(Neu2,2))*
     &           (SA*ZNeuC(Neu2,3) + CA*ZNeuC(Neu2,4))))/(CW*SW)

	ENDLOOP(Neu2)

#ifdef DETAILED_DEBUG
	DHIGGS "td_h0_mfv =", td_h0_mfv ENDL
#endif

	end


************************************************************************


	ComplexType function td_h0_nmfv(bX)
	implicit none
	integer bX

#include "FH.h"
#include "looptools.h"

	integer All2, Cha2, Gen2, Ind1, Ind2, Neu2

	ComplexType tmp1, tmp2, tmp3

	td_h0_nmfv = 0

	LOOP(Gen2, gM,3,1)

        td_h0_nmfv = td_h0_nmfv - 
     &    3/(8.D0*Pi**2)*(CA*EL1L*A0(Mf2(tM1,Gen2))*Mf2(tM1,Gen2))/
     &      (MW*SB*SW)

	ENDLOOP(Gen2)

	LOOP(All2, 1,6,gM)

	tmp1 = A0(MASf2(All2,tM1))

	LOOP(Ind2, 1,3,1)
	LOOP(Ind1, 1,3,1)

        td_h0_nmfv = td_h0_nmfv - 
     &    1/(32.D0*Pi**2)*(EL1L*tmp1*
     &        ((Delta(Ind1,Ind2)*
     &              (MW*MZ*SAB*SB*(3 - 4*SW2) - 
     &                6*CA*CW*Mf2(tM1,Ind1))*UASf(All2,Ind2,tM1) - 
     &             3*CW*(CA*Kf(Ind1,Ind2,tM1) + 
     &                MUEC*SA*Delta(Ind1,Ind2)*Mf(tM1,Ind1))*
     &              UASf(All2,3 + Ind2,tM1))*UASfC(All2,Ind1,tM1)-
     &            (3*CW*(CA*KfC(Ind2,Ind1,tM1) + 
     &                MUE*SA*Delta(Ind1,Ind2)*Mf(tM1,Ind1))*
     &              UASf(All2,Ind2,tM1) - 
     &             2*Delta(Ind1,Ind2)*
     &              (2*MW*MZ*SAB*SB*SW2 - 3*CA*CW*Mf2(tM1,Ind1))*
     &              UASf(All2,3 + Ind2,tM1))*
     &           UASfC(All2,3 + Ind1,tM1)))/(CW*MW*SB*SW)

	ENDLOOP(Ind1)
	ENDLOOP(Ind2)

	ENDLOOP(All2)

#ifdef DETAILED_DEBUG
	DHIGGS "td_h0_nmfv =", td_h0_nmfv ENDL
#endif

	if( mssmpart .le. 1 ) return

	LOOP(Gen2, gM,3,1)

        td_h0_nmfv = td_h0_nmfv + 
     &    3/(8.D0*Pi**2)*(EL1L*SA*A0(Mf2(bX,Gen2))*Mf2(bX,Gen2))/
     &      (CB*MW*SW)

	ENDLOOP(Gen2)

	LOOP(All2, 1,6,gM)

	tmp2 = A0(MASf2(All2,bX))

	LOOP(Ind2, 1,3,1)
	LOOP(Ind1, 1,3,1)

        td_h0_nmfv = td_h0_nmfv + 
     &    1/(32.D0*Pi**2)*(EL1L*tmp2*
     &        ((Delta(Ind1,Ind2)*
     &              (CB*MW*MZ*SAB*(3 - 2*SW2) - 
     &                6*CW*SA*Mf2(bX,Ind1))*UASf(All2,Ind2,bX) - 
     &             3*CW*(SA*Kf(Ind1,Ind2,bX) + 
     &                CA*MUEC*Delta(Ind1,Ind2)*Mf(bX,Ind1))*
     &              UASf(All2,3 + Ind2,bX))*UASfC(All2,Ind1,bX) - 
     &          (3*CW*(SA*KfC(Ind2,Ind1,bX) + 
     &                CA*MUE*Delta(Ind1,Ind2)*Mf(bX,Ind1))*
     &              UASf(All2,Ind2,bX) - 
     &             2*Delta(Ind1,Ind2)*
     &              (CB*MW*MZ*SAB*SW2 - 3*CW*SA*Mf2(bX,Ind1))*
     &              UASf(All2,3 + Ind2,bX))*UASfC(All2,3 + Ind1,bX)
     &          ))/(CB*CW*MW*SW)

	ENDLOOP(Ind1)
	ENDLOOP(Ind2)

	ENDLOOP(All2)

#ifdef DETAILED_DEBUG
	DHIGGS "td_h0_nmfv =", td_h0_nmfv ENDL
#endif

	if( mssmpart .le. 2 ) return

	LOOP(Gen2, 1,3,1)

        td_h0_nmfv = td_h0_nmfv - 
     &    1/(32.D0*Pi**2)*(EL1L*
     &        (CB*MW*MZ*SAB*A0(MASf2(Gen2,1)) - 
     &          4*CW*SA*A0(Mf2(2,Gen2))*Mf2(2,Gen2)))/(CB*CW*MW*SW)

	ENDLOOP(Gen2)

	LOOP(All2, 1,6,1)

	tmp3 = A0(MASf2(All2,2))

	LOOP(Ind2, 1,3,1)
	LOOP(Ind1, 1,3,1)

        td_h0_nmfv = td_h0_nmfv + 
     &    1/(32.D0*Pi**2)*(EL1L*tmp3*
     &        ((Delta(Ind1,Ind2)*
     &              (CB*MW*MZ*SAB*(1 - 2*SW2) - 
     &                2*CW*SA*Mf2(2,Ind1))*UASf(All2,Ind2,2) - 
     &             CW*(SA*Kf(Ind1,Ind2,2) + 
     &                CA*MUEC*Delta(Ind1,Ind2)*Mf(2,Ind1))*
     &              UASf(All2,3 + Ind2,2))*UASfC(All2,Ind1,2) - 
     &          (CW*(SA*KfC(Ind2,Ind1,2) + 
     &                CA*MUE*Delta(Ind1,Ind2)*Mf(2,Ind1))*
     &              UASf(All2,Ind2,2) - 
     &             2*Delta(Ind1,Ind2)*
     &              (CB*MW*MZ*SAB*SW2 - CW*SA*Mf2(2,Ind1))*
     &              UASf(All2,3 + Ind2,2))*UASfC(All2,3 + Ind1,2)))
     &       /(CB*CW*MW*SW)

	ENDLOOP(Ind1)
	ENDLOOP(Ind2)

	ENDLOOP(All2)

#ifdef DETAILED_DEBUG
	DHIGGS "td_h0_nmfv =", td_h0_nmfv ENDL
#endif

	if( mssmpart .le. 3 ) return

        td_h0_nmfv = td_h0_nmfv + 
     &    1/(64.D0*Pi**2)*(EL1L*MW*
     &        (SAB*(C2B*A0(MA02) + 3*C2A*A0(Mh02)) - 
     &          (2*CAB*S2A + C2A*SAB)*A0(MHH2) + 
     &          2*(C2B*SAB + 2*CW2*SBA)*A0(MHp2) - 
     &          2*(C2B*SAB - 6*CW2*SBA)*A0(MW2) - 
     &          (C2B*SAB - 6*SBA)*A0(MZ2)))/(CW2*SW)

	LOOP(Cha2, 1,2,1)

        td_h0_nmfv = td_h0_nmfv + 
     &    1/(8.D0*Pi**2*sqrt2)*
     &     (EL1L*A0(MCha2(Cha2))*MCha(Cha2)*
     &        (SA*(UCha(Cha2,2)*VCha(Cha2,1) + 
     &             UChaC(Cha2,2)*VChaC(Cha2,1)) - 
     &          CA*(UCha(Cha2,1)*VCha(Cha2,2) + 
     &             UChaC(Cha2,1)*VChaC(Cha2,2))))/SW

	ENDLOOP(Cha2)

	LOOP(Neu2, 1,4,1)

        td_h0_nmfv = td_h0_nmfv - 
     &    1/(16.D0*Pi**2)*(EL1L*A0(MNeu2(Neu2))*MNeu(Neu2)*
     &        ((SW*ZNeu(Neu2,1) - CW*ZNeu(Neu2,2))*
     &           (SA*ZNeu(Neu2,3) + CA*ZNeu(Neu2,4)) + 
     &          (SW*ZNeuC(Neu2,1) - CW*ZNeuC(Neu2,2))*
     &           (SA*ZNeuC(Neu2,3) + CA*ZNeuC(Neu2,4))))/(CW*SW)

	ENDLOOP(Neu2)

#ifdef DETAILED_DEBUG
	DHIGGS "td_h0_nmfv =", td_h0_nmfv ENDL
#endif

	end

