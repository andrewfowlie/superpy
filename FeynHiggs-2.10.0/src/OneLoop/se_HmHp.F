* se_HmHp.F
* this file is part of FeynHiggs
* generated 05 Nov 2013 19:22

#include "externals.h"
#include "types.h"
#include "debug.h"


	ComplexType function se_HmHp(k2, bX)
	implicit none
	RealType k2
	integer bX


#include "FH.h"

	ComplexType se_HmHp_mfv, se_HmHp_nmfv
	external se_HmHp_mfv, se_HmHp_nmfv

	if( ibits(fv, 3, 2) .eq. 0 ) then
	  se_HmHp = se_HmHp_mfv(k2, bX)
	else
	  se_HmHp = se_HmHp_nmfv(k2, bX)
	endif

#ifdef DEBUG
	DHIGGS 'k2 =', k2 ENDL
	DHIGGS 'se_HmHp =', se_HmHp ENDL
#endif
	end


************************************************************************


	ComplexType function se_HmHp_mfv(k2,bX)
	implicit none
	integer bX
	RealType k2

#include "FH.h"
#include "looptools.h"

	integer Cha3, Gen3, Gen4, Ind1, Neu3, Sfe3, Sfe4

	ComplexType dup1, dup2, dup3, dup4, dup5, tmp1, tmp2, tmp3
	ComplexType tmp4, tmp5, tmp6, tmp7

	se_HmHp_mfv = 0

	LOOP(Gen3, gM,3,1)
	LOOP(Sfe3, 1,2,1)

        se_HmHp_mfv = se_HmHp_mfv + 
     &    1/(16.D0*Pi)*(Alfa1L*
     &        (C2B*A0(MSf2(Sfe3,bX,Gen3))*
     &           (USf2(Sfe3,1,bX,Gen3) - 
     &             4*CW2*USf2(Sfe3,1,bX,Gen3) + 
     &             2*SW2*USf2(Sfe3,2,bX,Gen3)) + 
     &          A0(MSf2(Sfe3,tM1,Gen3))*
     &           (C2B*(1 + 2*CW2)*USf2(Sfe3,1,tM1,Gen3) - 
     &             (2*(2*C2B*MW2*SW2*TB2 - 3*CW2*Mf2(tM1,Gen3))*
     &                USf2(Sfe3,2,tM1,Gen3))/(MW2*TB2))))/(CW2*SW2)

	ENDLOOP(Sfe3)
	ENDLOOP(Gen3)

	LOOP(Gen3, gM,3,1)

	tmp1 = B0(k2,0.D0,Mf2(tM1,Gen3))

	tmp2 = B1(k2,0.D0,Mf2(tM1,Gen3))

	LOOP(Gen4, gM,3,1)

        se_HmHp_mfv = se_HmHp_mfv + 
     &    3/(4.D0*Pi)*(Alfa1L*CKMin(Gen3,Gen4)*CKMinC(Gen3,Gen4)*
     &        (k2*tmp2 + tmp1*(k2 - Mf2(tM1,Gen3)))*Mf2(tM1,Gen3))/
     &      (MW2*SW2*TB2)

	ENDLOOP(Gen4)

	ENDLOOP(Gen3)

	LOOP(Gen3, gM,3,1)
	LOOP(Sfe3, 1,2,1)

	tmp3 = A0(MSf2(Sfe3,bX,Gen3))

	LOOP(Ind1, 1,3,1)

        se_HmHp_mfv = se_HmHp_mfv + 
     &    3/(8.D0*Pi)*(Alfa1L*tmp3*CKMin(Ind1,Gen3)*
     &        CKMinC(Ind1,Gen3)*Mf2(tM1,Ind1)*USf2(Sfe3,1,bX,Gen3))
     &       /(MW2*SW2*TB2)

	ENDLOOP(Ind1)

	ENDLOOP(Sfe3)
	ENDLOOP(Gen3)

	LOOP(Gen4, gM,3,1)
	LOOP(Gen3, gM,3,1)
	LOOP(Sfe4, 1,2,1)
	LOOP(Sfe3, 1,2,1)

        se_HmHp_mfv = se_HmHp_mfv + 
     &    3/(8.D0*Pi)*(Alfa1L*
     &        B0(k2,MSf2(Sfe3,tM1,Gen3),MSf2(Sfe4,bX,Gen4))*
     &        CKMin(Gen3,Gen4)*CKMinC(Gen3,Gen4)*
     &        ((MW2*S2B*TB - Mf2(tM1,Gen3))*USf(Sfe3,1,tM1,Gen3) - 
     &          (Kf(Gen3,Gen3,tM1) + MUEC*TB*Mf(tM1,Gen3))*
     &           USf(Sfe3,2,tM1,Gen3))*USf(Sfe4,1,bX,Gen4)*
     &        ((MW2*S2B*TB - Mf2(tM1,Gen3))*
     &           USfC(Sfe3,1,tM1,Gen3) - 
     &          (KfC(Gen3,Gen3,tM1) + MUE*TB*Mf(tM1,Gen3))*
     &           USfC(Sfe3,2,tM1,Gen3))*USfC(Sfe4,1,bX,Gen4))/
     &      (MW2*SW2*TB2)

	ENDLOOP(Sfe3)
	ENDLOOP(Sfe4)
	ENDLOOP(Gen3)
	ENDLOOP(Gen4)

#ifdef DETAILED_DEBUG
	DHIGGS "se_HmHp_mfv =", se_HmHp_mfv ENDL
#endif

	if( mssmpart .le. 1 ) return

	LOOP(Gen3, gM,3,1)
	LOOP(Sfe3, 1,2,1)

        se_HmHp_mfv = se_HmHp_mfv + 
     &    3/(8.D0*Pi)*(Alfa1L*TB2*A0(MSf2(Sfe3,bX,Gen3))*
     &        Mf2(bX,Gen3)*USf2(Sfe3,2,bX,Gen3))/(MW2*SW2)

	ENDLOOP(Sfe3)
	ENDLOOP(Gen3)

	LOOP(Gen3, gM,3,1)

	tmp4 = B0(k2,0.D0,Mf2(tM1,Gen3))

	tmp5 = B1(k2,0.D0,Mf2(tM1,Gen3))

	LOOP(Gen4, gM,3,1)

        se_HmHp_mfv = se_HmHp_mfv - 
     &    3/(4.D0*Pi)*(Alfa1L*CKMin(Gen3,Gen4)*CKMinC(Gen3,Gen4)*
     &        ((A0(Mf2(bX,Gen4)) - 
     &             k2*B1(k2,Mf2(bX,Gen4),Mf2(tM1,Gen3)))*
     &           (TB2**2*Mf2(bX,Gen4) + Mf2(tM1,Gen3)) + 
     &          Mf2(tM1,Gen3)*
     &           (k2*(tmp4 + tmp5) - tmp4*Mf2(tM1,Gen3)) - 
     &          B0(k2,Mf2(bX,Gen4),Mf2(tM1,Gen3))*
     &           ((k2 - Mf2(tM1,Gen3))*Mf2(tM1,Gen3) + 
     &             TB2*Mf2(bX,Gen4)*
     &              (k2*TB2 + (-2 - TB2)*Mf2(tM1,Gen3)))))/
     &      (MW2*SW2*TB2)

	ENDLOOP(Gen4)

	ENDLOOP(Gen3)

	LOOP(Gen3, gM,3,1)
	LOOP(Sfe3, 1,2,1)

	tmp6 = A0(MSf2(Sfe3,tM1,Gen3))

	LOOP(Ind1, 1,3,1)

        se_HmHp_mfv = se_HmHp_mfv + 
     &    3/(8.D0*Pi)*(Alfa1L*TB2*tmp6*CKMin(Gen3,Ind1)*
     &        CKMinC(Gen3,Ind1)*Mf2(bX,Ind1)*USf2(Sfe3,1,tM1,Gen3))
     &       /(MW2*SW2)

	ENDLOOP(Ind1)

	ENDLOOP(Sfe3)
	ENDLOOP(Gen3)

	LOOP(Gen4, gM,3,1)
	LOOP(Gen3, gM,3,1)
	LOOP(Sfe4, 1,2,1)
	LOOP(Sfe3, 1,2,1)

	dup1 = KfC(Gen3,Gen3,tM1) + MUE*TB*Mf(tM1,Gen3)

	dup2 = Kf(Gen3,Gen3,tM1) + MUEC*TB*Mf(tM1,Gen3)

	dup3 = -(MW2*S2B*TB) + TB2*Mf2(bX,Gen4) + Mf2(tM1,Gen3)

        se_HmHp_mfv = se_HmHp_mfv - 
     &    3/(8.D0*Pi)*(Alfa1L*
     &        B0(k2,MSf2(Sfe3,tM1,Gen3),MSf2(Sfe4,bX,Gen4))*
     &        CKMin(Gen3,Gen4)*CKMinC(Gen3,Gen4)*
     &        (((MW2*S2B*TB - Mf2(tM1,Gen3))*
     &              USf(Sfe3,1,tM1,Gen3) - 
     &             dup2*USf(Sfe3,2,tM1,Gen3))*USf(Sfe4,1,bX,Gen4)*
     &           ((MW2*S2B*TB - Mf2(tM1,Gen3))*
     &              USfC(Sfe3,1,tM1,Gen3) - 
     &             dup1*USfC(Sfe3,2,tM1,Gen3))*USfC(Sfe4,1,bX,Gen4)
     &            - (USf(Sfe4,1,bX,Gen4)*
     &              (dup3*USfC(Sfe3,1,tM1,Gen3) + 
     &                dup1*USfC(Sfe3,2,tM1,Gen3)) + 
     &             USf(Sfe4,2,bX,Gen4)*
     &              ((TB2*Kf(Gen4,Gen4,bX) + MUEC*TB*Mf(bX,Gen4))*
     &                 USfC(Sfe3,1,tM1,Gen3) + 
     &                (1 + TB2)*Mf(bX,Gen4)*Mf(tM1,Gen3)*
     &                 USfC(Sfe3,2,tM1,Gen3)))*
     &           (USf(Sfe3,1,tM1,Gen3)*
     &              (dup3*USfC(Sfe4,1,bX,Gen4) + 
     &                (TB2*KfC(Gen4,Gen4,bX) + MUE*TB*Mf(bX,Gen4))*
     &                 USfC(Sfe4,2,bX,Gen4)) + 
     &             USf(Sfe3,2,tM1,Gen3)*
     &              (dup2*USfC(Sfe4,1,bX,Gen4) + 
     &                Mf(bX,Gen4)*Mf(tM1,Gen3)*
     &                 (USfC(Sfe4,2,bX,Gen4) + 
     &                   TB2*USfC(Sfe4,2,bX,Gen4))))))/
     &      (MW2*SW2*TB2)

	ENDLOOP(Sfe3)
	ENDLOOP(Sfe4)
	ENDLOOP(Gen3)
	ENDLOOP(Gen4)

#ifdef DETAILED_DEBUG
	DHIGGS "se_HmHp_mfv =", se_HmHp_mfv ENDL
#endif

	if( mssmpart .le. 2 ) return

	LOOP(Gen3, 1,3,1)

        se_HmHp_mfv = se_HmHp_mfv - 
     &    1/(16.D0*Pi)*(Alfa1L*
     &        (4*CW2*TB2*(A0(Mf2(2,Gen3)) + 
     &             k2*B1(k2,0.D0,Mf2(2,Gen3)))*Mf2(2,Gen3) + 
     &          A0(MSf2(1,1,Gen3))*
     &           (C2B*(MW2 - 2*CW2*MW2) - 2*CW2*TB2*Mf2(2,Gen3))))/
     &      (CW2*MW2*SW2)

	ENDLOOP(Gen3)

	LOOP(Sfe3, 1,2,1)
	LOOP(Gen3, 1,3,1)

        se_HmHp_mfv = se_HmHp_mfv - 
     &    1/(16.D0*Pi)*(Alfa1L*
     &        (A0(MSf2(Sfe3,2,Gen3))*
     &           (C2B*MW2*USf2(Sfe3,1,2,Gen3) - 
     &             2*(C2B*MW2*SW2 + CW2*TB2*Mf2(2,Gen3))*
     &              USf2(Sfe3,2,2,Gen3)) - 
     &          2*CW2*B0(k2,MSf2(1,1,Gen3),MSf2(Sfe3,2,Gen3))*
     &           ((MW2*S2B - TB*Mf2(2,Gen3))*USf(Sfe3,1,2,Gen3) - 
     &             (TB*Kf(Gen3,Gen3,2) + MUEC*Mf(2,Gen3))*
     &              USf(Sfe3,2,2,Gen3))*
     &           ((MW2*S2B - TB*Mf2(2,Gen3))*USfC(Sfe3,1,2,Gen3) - 
     &             (TB*KfC(Gen3,Gen3,2) + MUE*Mf(2,Gen3))*
     &              USfC(Sfe3,2,2,Gen3))))/(CW2*MW2*SW2)

	ENDLOOP(Gen3)
	ENDLOOP(Sfe3)

#ifdef DETAILED_DEBUG
	DHIGGS "se_HmHp_mfv =", se_HmHp_mfv ENDL
#endif

	if( mssmpart .le. 3 ) return

        se_HmHp_mfv = se_HmHp_mfv - 
     &    1/(32.D0*Pi)*(Alfa1L*
     &        (2*((CBA2*CW2**2*(k2 + Mh02 - MW2) + 
     &                MW2*(2*CBA*CW2*S2B*SAB - S2B**2*SAB**2))*
     &              B0(k2,Mh02,MW2) - 
     &             MW2*((C2B*(C2B*SAB**2 + 4*CW2*SAB*SBA) + 
     &                   4*CW2**2*SBA2)*B0(k2,Mh02,MHp2) + 
     &                (4*CBA2*CW2**2 + 
     &                   C2B*(C2B*CAB**2 - 4*CAB*CBA*CW2))*
     &                 B0(k2,MHH2,MHp2)) - 
     &             (CAB*MW2*(CAB*S2B**2 - 2*CW2*S2B*SBA) - 
     &                CW2**2*(k2 + MHH2 - MW2)*SBA2)*
     &              B0(k2,MHH2,MW2)) - 
     &          CW2*((CW2*(1 - S2A*S2B) + C2A*C2B*SW2)*A0(Mh02) + 
     &             (CW2*(1 + S2A*S2B) - C2A*C2B*SW2)*A0(MHH2) + 
     &             C2B**2*(A0(MA02) + 4*A0(MHp2)) - 
     &             2*(C2B**2 - 6*CW2 - S2B**2)*A0(MW2) + 
     &             (CW2*(1 + 6*CW2 + S2B**2 - 12*SW2) - 
     &                C2B**2*SW2 + 6*SW2**2)*A0(MZ2) - 
     &             CW2*(2*(k2 + MA02 - MW2)*B0(k2,MA02,MW2) - 
     &                4*k2*
     &                 (B1(k2,MA02,MW2) + CBA2*B1(k2,Mh02,MW2) + 
     &                   SBA2*B1(k2,MHH2,MW2)) + 
     &                SW2*(8*MHp2*B0(k2,0.D0,MHp2) + 
     &                   k2*(8*B0(k2,0.D0,MHp2) - 16*B1(k2,MHp2,0.D0)))
     &                ) - (CW2 - SW2)**2*
     &              (2*MHp2*B0(k2,MHp2,MZ2) + 
     &                k2*(2*B0(k2,MHp2,MZ2) - 4*B1(k2,MHp2,MZ2)))))
     &        )/(CW2**2*SW2)

	LOOP(Cha3, 1,2,1)

	tmp7 = A0(MCha2(Cha3))

	LOOP(Neu3, 1,4,1)

	dup4 = (SW*ZNeu(Neu3,1) + CW*ZNeu(Neu3,2))/CW

	dup5 = (SW*ZNeuC(Neu3,1) + CW*ZNeuC(Neu3,2))/CW

        se_HmHp_mfv = se_HmHp_mfv - 
     &    1/(8.D0*Pi)*(Alfa1L*
     &        ((tmp7 - k2*B1(k2,MCha2(Cha3),MNeu2(Neu3)))*
     &           (4/sqrt2**2*
     &              (CB2*(sqrt2*(VCha(Cha3,1)*ZNeu(Neu3,4)) + 
     &                  dup4*VCha(Cha3,2))*
     &                (sqrt2*(VChaC(Cha3,1)*ZNeuC(Neu3,4)) + 
     &                  dup5*VChaC(Cha3,2))) + 
     &             SB2*(sqrt2*(dup4*UCha(Cha3,2)) - 
     &                2*UCha(Cha3,1)*ZNeu(Neu3,3))*
     &              (sqrt2*(dup5*UChaC(Cha3,2)) - 
     &                2*UChaC(Cha3,1)*ZNeuC(Neu3,3))) - 
     &          B0(k2,MCha2(Cha3),MNeu2(Neu3))*
     &           (4/sqrt2**2*
     &              (CB*SB*MCha(Cha3)*MNeu(Neu3)*
     &                ((-(sqrt2*(UCha(Cha3,1)*ZNeu(Neu3,3))) + 
     &                     dup4*UCha(Cha3,2))*
     &                   (sqrt2*(VCha(Cha3,1)*ZNeu(Neu3,4)) + 
     &                     dup4*VCha(Cha3,2)) + 
     &                  (-(sqrt2*(UChaC(Cha3,1)*ZNeuC(Neu3,3))) + 
     &                     dup5*UChaC(Cha3,2))*
     &                   (sqrt2*(VChaC(Cha3,1)*ZNeuC(Neu3,4)) + 
     &                     dup5*VChaC(Cha3,2)))) + 
     &             (k2 - MNeu2(Neu3))*
     &              (4/sqrt2**2*
     &                 (CB2*
     &                   (sqrt2*(VCha(Cha3,1)*ZNeu(Neu3,4)) + 
     &                     dup4*VCha(Cha3,2))*
     &                   (sqrt2*(VChaC(Cha3,1)*ZNeuC(Neu3,4)) + 
     &                     dup5*VChaC(Cha3,2))) + 
     &                SB2*(sqrt2*(dup4*UCha(Cha3,2)) - 
     &                   2*UCha(Cha3,1)*ZNeu(Neu3,3))*
     &                 (sqrt2*(dup5*UChaC(Cha3,2)) - 
     &                   2*UChaC(Cha3,1)*ZNeuC(Neu3,3))))))/SW2

	ENDLOOP(Neu3)

	ENDLOOP(Cha3)

#ifdef DETAILED_DEBUG
	DHIGGS "se_HmHp_mfv =", se_HmHp_mfv ENDL
#endif

	end


************************************************************************


	ComplexType function se_HmHp_nmfv(k2,bX)
	implicit none
	integer bX
	RealType k2

#include "FH.h"
#include "looptools.h"

	integer All3, All4, Cha3, Gen3, Gen4, Ind1, Ind2, Ind3, Ind4
	integer Ind5, Ind6, Neu3

	ComplexType dup1, dup2, dup3, tmp1, tmp10, tmp11, tmp12
	ComplexType tmp13, tmp14, tmp15, tmp16, tmp17, tmp18, tmp19
	ComplexType tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9

	se_HmHp_nmfv = 0

	LOOP(All3, 1,6,gM)

	tmp1 = A0(MASf2(All3,bX))

	tmp2 = A0(MASf2(All3,tM1))

	LOOP(Ind1, 1,3,1)

        se_HmHp_nmfv = se_HmHp_nmfv + 
     &    1/(16.D0*Pi)*(Alfa1L*
     &        (C2B*tmp1*(UASf(All3,Ind1,bX)*
     &              (UASfC(All3,Ind1,bX) - 
     &                4*CW2*UASfC(All3,Ind1,bX)) + 
     &             2*SW2*UASf(All3,3 + Ind1,bX)*
     &              UASfC(All3,3 + Ind1,bX)) + 
     &          tmp2*(C2B*(1 + 2*CW2)*UASf(All3,Ind1,tM1)*
     &              UASfC(All3,Ind1,tM1) - 
     &             (2*(2*C2B*MW2*SW2*TB2 - 3*CW2*Mf2(tM1,Ind1))*
     &                UASf(All3,3 + Ind1,tM1)*
     &                UASfC(All3,3 + Ind1,tM1))/(MW2*TB2))))/
     &      (CW2*SW2)

	ENDLOOP(Ind1)

	ENDLOOP(All3)

	LOOP(Gen3, gM,3,1)

	tmp3 = B0(k2,0.D0,Mf2(tM1,Gen3))

	tmp4 = B1(k2,0.D0,Mf2(tM1,Gen3))

	LOOP(Gen4, gM,3,1)

        se_HmHp_nmfv = se_HmHp_nmfv + 
     &    3/(4.D0*Pi)*(Alfa1L*CKMin(Gen3,Gen4)*CKMinC(Gen3,Gen4)*
     &        (k2*tmp4 + tmp3*(k2 - Mf2(tM1,Gen3)))*Mf2(tM1,Gen3))/
     &      (MW2*SW2*TB2)

	ENDLOOP(Gen4)

	ENDLOOP(Gen3)

	LOOP(All3, 1,6,gM)

	tmp5 = A0(MASf2(All3,bX))

	LOOP(Ind3, 1,3,1)
	LOOP(Ind2, 1,3,1)
	LOOP(Ind1, 1,3,1)

        se_HmHp_nmfv = se_HmHp_nmfv + 
     &    3/(8.D0*Pi)*(Alfa1L*tmp5*CKMin(Ind1,Ind2)*
     &        CKMinC(Ind1,Ind3)*Mf2(tM1,Ind1)*UASf(All3,Ind3,bX)*
     &        UASfC(All3,Ind2,bX))/(MW2*SW2*TB2)

	ENDLOOP(Ind1)
	ENDLOOP(Ind2)
	ENDLOOP(Ind3)

	ENDLOOP(All3)

	LOOP(All4, 1,6,gM)
	LOOP(All3, 1,6,gM)

	tmp6 = B0(k2,MASf2(All3,tM1),MASf2(All4,bX))

	LOOP(Ind4, 1,3,1)
	LOOP(Ind3, 1,3,1)
	LOOP(Ind2, 1,3,1)
	LOOP(Ind1, 1,3,1)

        se_HmHp_nmfv = se_HmHp_nmfv - 
     &    3/(8.D0*Pi)*(Alfa1L*tmp6*CKMin(Ind1,Ind2)*
     &        CKMinC(Ind3,Ind4)*UASf(All4,Ind4,bX)*
     &        (MUEC*Mf(tM1,Ind1)*UASf(All3,3 + Ind1,tM1)*
     &           ((MW2*S2B*TB2 - TB*Mf2(tM1,Ind3))*
     &              UASfC(All3,Ind3,tM1) - 
     &             MUE*TB2*Mf(tM1,Ind3)*UASfC(All3,3 + Ind3,tM1))-
     &            UASf(All3,Ind1,tM1)*
     &           ((MW2*S2B*(MW2*S2B*TB2 - TB*Mf2(tM1,Ind1)) - 
     &                (MW2*S2B*TB - Mf2(tM1,Ind1))*Mf2(tM1,Ind3))*
     &              UASfC(All3,Ind3,tM1) - 
     &             MUE*Mf(tM1,Ind3)*
     &              (MW2*S2B*TB2 - TB*Mf2(tM1,Ind1))*
     &              UASfC(All3,3 + Ind3,tM1)))*UASfC(All4,Ind2,bX))
     &       /(MW2*SW2*TB2)

	ENDLOOP(Ind1)
	ENDLOOP(Ind2)
	ENDLOOP(Ind3)
	ENDLOOP(Ind4)

	ENDLOOP(All3)
	ENDLOOP(All4)

	LOOP(All4, 1,6,gM)
	LOOP(All3, 1,6,gM)

	tmp7 = B0(k2,MASf2(All3,tM1),MASf2(All4,bX))

	LOOP(Ind5, 1,3,1)
	LOOP(Ind4, 1,3,1)
	LOOP(Ind3, 1,3,1)
	LOOP(Ind2, 1,3,1)
	LOOP(Ind1, 1,3,1)

        se_HmHp_nmfv = se_HmHp_nmfv - 
     &    3/(8.D0*Pi)*(Alfa1L*tmp7*CKMin(Ind1,Ind2)*
     &        CKMinC(Ind3,Ind4)*UASf(All4,Ind4,bX)*
     &        (Kf(Ind1,Ind5,tM1)*UASf(All3,3 + Ind5,tM1)*
     &           ((MW2*S2B*TB - Mf2(tM1,Ind3))*
     &              UASfC(All3,Ind3,tM1) - 
     &             MUE*TB*Mf(tM1,Ind3)*UASfC(All3,3 + Ind3,tM1)) + 
     &          KfC(Ind3,Ind5,tM1)*
     &           ((MW2*S2B*TB - Mf2(tM1,Ind1))*
     &              UASf(All3,Ind1,tM1) - 
     &             MUEC*TB*Mf(tM1,Ind1)*UASf(All3,3 + Ind1,tM1))*
     &           UASfC(All3,3 + Ind5,tM1))*UASfC(All4,Ind2,bX))/
     &      (MW2*SW2*TB2)

	ENDLOOP(Ind1)
	ENDLOOP(Ind2)
	ENDLOOP(Ind3)
	ENDLOOP(Ind4)
	ENDLOOP(Ind5)

	ENDLOOP(All3)
	ENDLOOP(All4)

	LOOP(All4, 1,6,gM)
	LOOP(All3, 1,6,gM)

	tmp8 = B0(k2,MASf2(All3,tM1),MASf2(All4,bX))

	LOOP(Ind6, 1,3,1)
	LOOP(Ind5, 1,3,1)
	LOOP(Ind4, 1,3,1)
	LOOP(Ind3, 1,3,1)
	LOOP(Ind2, 1,3,1)
	LOOP(Ind1, 1,3,1)

        se_HmHp_nmfv = se_HmHp_nmfv + 
     &    3/(8.D0*Pi)*(Alfa1L*tmp8*CKMin(Ind1,Ind2)*
     &        CKMinC(Ind3,Ind4)*Kf(Ind1,Ind5,tM1)*
     &        KfC(Ind3,Ind6,tM1)*UASf(All3,3 + Ind5,tM1)*
     &        UASf(All4,Ind4,bX)*UASfC(All3,3 + Ind6,tM1)*
     &        UASfC(All4,Ind2,bX))/(MW2*SW2*TB2)

	ENDLOOP(Ind1)
	ENDLOOP(Ind2)
	ENDLOOP(Ind3)
	ENDLOOP(Ind4)
	ENDLOOP(Ind5)
	ENDLOOP(Ind6)

	ENDLOOP(All3)
	ENDLOOP(All4)

#ifdef DETAILED_DEBUG
	DHIGGS "se_HmHp_nmfv =", se_HmHp_nmfv ENDL
#endif

	if( mssmpart .le. 1 ) return

	LOOP(All3, 1,6,gM)

	tmp9 = A0(MASf2(All3,bX))

	LOOP(Ind1, 1,3,1)

        se_HmHp_nmfv = se_HmHp_nmfv + 
     &    3/(8.D0*Pi)*(Alfa1L*TB2*tmp9*Mf2(bX,Ind1)*
     &        UASf(All3,3 + Ind1,bX)*UASfC(All3,3 + Ind1,bX))/
     &      (MW2*SW2)

	ENDLOOP(Ind1)

	ENDLOOP(All3)

	LOOP(Gen3, gM,3,1)

	tmp10 = B0(k2,0.D0,Mf2(tM1,Gen3))

	tmp11 = B1(k2,0.D0,Mf2(tM1,Gen3))

	LOOP(Gen4, gM,3,1)

        se_HmHp_nmfv = se_HmHp_nmfv - 
     &    3/(4.D0*Pi)*(Alfa1L*CKMin(Gen3,Gen4)*CKMinC(Gen3,Gen4)*
     &        ((A0(Mf2(bX,Gen4)) - 
     &             k2*B1(k2,Mf2(bX,Gen4),Mf2(tM1,Gen3)))*
     &           (TB2**2*Mf2(bX,Gen4) + Mf2(tM1,Gen3)) + 
     &          Mf2(tM1,Gen3)*
     &           (k2*(tmp10 + tmp11) - tmp10*Mf2(tM1,Gen3)) - 
     &          B0(k2,Mf2(bX,Gen4),Mf2(tM1,Gen3))*
     &           ((k2 - Mf2(tM1,Gen3))*Mf2(tM1,Gen3) + 
     &             TB2*Mf2(bX,Gen4)*
     &              (k2*TB2 + (-2 - TB2)*Mf2(tM1,Gen3)))))/
     &      (MW2*SW2*TB2)

	ENDLOOP(Gen4)

	ENDLOOP(Gen3)

	LOOP(All3, 1,6,gM)

	tmp12 = A0(MASf2(All3,tM1))

	LOOP(Ind3, 1,3,1)
	LOOP(Ind2, 1,3,1)
	LOOP(Ind1, 1,3,1)

        se_HmHp_nmfv = se_HmHp_nmfv + 
     &    3/(8.D0*Pi)*(Alfa1L*TB2*tmp12*CKMin(Ind1,Ind2)*
     &        CKMinC(Ind3,Ind2)*Mf2(bX,Ind2)*UASf(All3,Ind1,tM1)*
     &        UASfC(All3,Ind3,tM1))/(MW2*SW2)

	ENDLOOP(Ind1)
	ENDLOOP(Ind2)
	ENDLOOP(Ind3)

	ENDLOOP(All3)

	LOOP(All4, 1,6,gM)
	LOOP(All3, 1,6,gM)

	tmp13 = B0(k2,MASf2(All3,tM1),MASf2(All4,bX))

	LOOP(Ind4, 1,3,1)
	LOOP(Ind3, 1,3,1)
	LOOP(Ind2, 1,3,1)
	LOOP(Ind1, 1,3,1)

        se_HmHp_nmfv = se_HmHp_nmfv - 
     &    3/(8.D0*Pi)*(Alfa1L*tmp13*CKMin(Ind1,Ind2)*
     &        CKMinC(Ind3,Ind4)*
     &        (((TB2*Mf2(bX,Ind2)*
     &                 (MW2*S2B*TB - TB2*Mf2(bX,Ind4) - 
     &                   Mf2(tM1,Ind3))*UASf(All3,Ind1,tM1)*
     &                 UASf(All4,Ind4,bX) + 
     &                MUEC*Mf(bX,Ind4)*
     &                 ((TB2*(MW2*S2B - TB*Mf2(bX,Ind2)) - 
     &                      TB*Mf2(tM1,Ind1))*UASf(All3,Ind1,tM1)-
     &                     MUEC*TB2*Mf(tM1,Ind1)*
     &                    UASf(All3,3 + Ind1,tM1))*
     &                 UASf(All4,3 + Ind4,bX))*UASfC(All3,Ind3,tM1)
     &               - Mf(tM1,Ind3)*
     &              (MUE*TB*TB2*Mf2(bX,Ind2)*UASf(All3,Ind1,tM1)*
     &                 UASf(All4,Ind4,bX) - 
     &                (1 + TB2)*Mf(bX,Ind4)*
     &                 ((MW2*S2B*TB - TB2*Mf2(bX,Ind2) - 
     &                      Mf2(tM1,Ind1))*UASf(All3,Ind1,tM1) - 
     &                   MUEC*TB*Mf(tM1,Ind1)*
     &                    UASf(All3,3 + Ind1,tM1))*
     &                 UASf(All4,3 + Ind4,bX))*
     &              UASfC(All3,3 + Ind3,tM1))*UASfC(All4,Ind2,bX)+
     &            Mf(bX,Ind2)*
     &           ((1 + TB2)*Mf(tM1,Ind1)*UASf(All3,3 + Ind1,tM1)*
     &              (((MW2*S2B*TB - Mf2(tM1,Ind3))*
     &                    UASf(All4,Ind4,bX) - 
     &                   MUEC*TB*Mf(bX,Ind4)*UASf(All4,3 + Ind4,bX)
     &                   )*UASfC(All3,Ind3,tM1) - 
     &                Mf(tM1,Ind3)*
     &                 (MUE*TB*UASf(All4,Ind4,bX) + 
     &                   Mf(bX,Ind4)*
     &                    (UASf(All4,3 + Ind4,bX) + 
     &                      TB2*UASf(All4,3 + Ind4,bX)))*
     &                 UASfC(All3,3 + Ind3,tM1)) + 
     &             MUE*UASf(All3,Ind1,tM1)*
     &              (((MW2*S2B*TB2 - TB*Mf2(tM1,Ind3))*
     &                    UASf(All4,Ind4,bX) - 
     &                   MUEC*TB2*Mf(bX,Ind4)*
     &                    UASf(All4,3 + Ind4,bX))*
     &                 UASfC(All3,Ind3,tM1) - 
     &                Mf(tM1,Ind3)*
     &                 (MUE*TB2*UASf(All4,Ind4,bX) + 
     &                   TB*Mf(bX,Ind4)*
     &                    (UASf(All4,3 + Ind4,bX) + 
     &                      TB2*UASf(All4,3 + Ind4,bX)))*
     &                 UASfC(All3,3 + Ind3,tM1)))*
     &           UASfC(All4,3 + Ind2,bX) + 
     &          TB2*Mf2(bX,Ind4)*UASf(All4,Ind4,bX)*
     &           UASfC(All3,Ind3,tM1)*
     &           (UASf(All3,Ind1,tM1)*
     &              ((MW2*S2B*TB - Mf2(tM1,Ind1))*
     &                 UASfC(All4,Ind2,bX) - 
     &                MUE*TB*Mf(bX,Ind2)*UASfC(All4,3 + Ind2,bX))-
     &               Mf(tM1,Ind1)*UASf(All3,3 + Ind1,tM1)*
     &              (MUEC*TB*UASfC(All4,Ind2,bX) + 
     &                Mf(bX,Ind2)*
     &                 (UASfC(All4,3 + Ind2,bX) + 
     &                   TB2*UASfC(All4,3 + Ind2,bX))))))/
     &      (MW2*SW2*TB2)

	ENDLOOP(Ind1)
	ENDLOOP(Ind2)
	ENDLOOP(Ind3)
	ENDLOOP(Ind4)

	ENDLOOP(All3)
	ENDLOOP(All4)

	LOOP(All4, 1,6,gM)
	LOOP(All3, 1,6,gM)

	tmp14 = B0(k2,MASf2(All3,tM1),MASf2(All4,bX))

	LOOP(Ind5, 1,3,1)
	LOOP(Ind4, 1,3,1)
	LOOP(Ind3, 1,3,1)
	LOOP(Ind2, 1,3,1)
	LOOP(Ind1, 1,3,1)

        dup1 = MUE*TB*UASf(All3,Ind1,tM1) + 
     &    Mf(tM1,Ind1)*(UASf(All3,3 + Ind1,tM1) + 
     &       TB2*UASf(All3,3 + Ind1,tM1))

        se_HmHp_nmfv = se_HmHp_nmfv + 
     &    3/(8.D0*Pi)*(Alfa1L*tmp14*CKMin(Ind1,Ind2)*
     &        CKMinC(Ind3,Ind4)*
     &        (Kf(Ind1,Ind5,tM1)*UASf(All3,3 + Ind5,tM1)*
     &           (TB2*Mf2(bX,Ind4)*UASf(All4,Ind4,bX)*
     &              UASfC(All3,Ind3,tM1) + 
     &             Mf(bX,Ind4)*UASf(All4,3 + Ind4,bX)*
     &              (MUEC*TB*UASfC(All3,Ind3,tM1) + 
     &                Mf(tM1,Ind3)*
     &                 (UASfC(All3,3 + Ind3,tM1) + 
     &                   TB2*UASfC(All3,3 + Ind3,tM1))))*
     &           UASfC(All4,Ind2,bX) + 
     &          KfC(Ind3,Ind5,tM1)*UASf(All4,Ind4,bX)*
     &           UASfC(All3,3 + Ind5,tM1)*
     &           (TB2*Mf2(bX,Ind2)*UASf(All3,Ind1,tM1)*
     &              UASfC(All4,Ind2,bX) + 
     &             dup1*Mf(bX,Ind2)*UASfC(All4,3 + Ind2,bX)) - 
     &          TB2*(Kf(Ind4,Ind5,bX)*UASf(All4,3 + Ind5,bX)*
     &              UASfC(All3,Ind3,tM1)*
     &              (((MW2*S2B*TB - TB2*Mf2(bX,Ind2) - 
     &                      Mf2(tM1,Ind1))*UASf(All3,Ind1,tM1) - 
     &                   MUEC*TB*Mf(tM1,Ind1)*
     &                    UASf(All3,3 + Ind1,tM1))*
     &                 UASfC(All4,Ind2,bX) - 
     &                dup1*Mf(bX,Ind2)*UASfC(All4,3 + Ind2,bX)) + 
     &             KfC(Ind2,Ind5,bX)*UASf(All3,Ind1,tM1)*
     &              (((MW2*S2B*TB - TB2*Mf2(bX,Ind4) - 
     &                      Mf2(tM1,Ind3))*UASf(All4,Ind4,bX) - 
     &                   MUEC*TB*Mf(bX,Ind4)*UASf(All4,3 + Ind4,bX)
     &                   )*UASfC(All3,Ind3,tM1) - 
     &                Mf(tM1,Ind3)*
     &                 (MUE*TB*UASf(All4,Ind4,bX) + 
     &                   Mf(bX,Ind4)*
     &                    (UASf(All4,3 + Ind4,bX) + 
     &                      TB2*UASf(All4,3 + Ind4,bX)))*
     &                 UASfC(All3,3 + Ind3,tM1))*
     &              UASfC(All4,3 + Ind5,bX))))/(MW2*SW2*TB2)

	ENDLOOP(Ind1)
	ENDLOOP(Ind2)
	ENDLOOP(Ind3)
	ENDLOOP(Ind4)
	ENDLOOP(Ind5)

	ENDLOOP(All3)
	ENDLOOP(All4)

	LOOP(All4, 1,6,gM)
	LOOP(All3, 1,6,gM)

	tmp15 = B0(k2,MASf2(All3,tM1),MASf2(All4,bX))

	LOOP(Ind6, 1,3,1)
	LOOP(Ind5, 1,3,1)
	LOOP(Ind4, 1,3,1)
	LOOP(Ind3, 1,3,1)
	LOOP(Ind2, 1,3,1)
	LOOP(Ind1, 1,3,1)

        se_HmHp_nmfv = se_HmHp_nmfv + 
     &    3/(8.D0*Pi)*(Alfa1L*tmp15*CKMin(Ind1,Ind2)*
     &        CKMinC(Ind3,Ind4)*
     &        (Kf(Ind1,Ind5,tM1)*Kf(Ind4,Ind6,bX)*
     &           UASf(All3,3 + Ind5,tM1)*UASf(All4,3 + Ind6,bX)*
     &           UASfC(All3,Ind3,tM1)*UASfC(All4,Ind2,bX) + 
     &          KfC(Ind2,Ind6,bX)*UASf(All3,Ind1,tM1)*
     &           (TB2*Kf(Ind4,Ind5,bX)*UASf(All4,3 + Ind5,bX)*
     &              UASfC(All3,Ind3,tM1) + 
     &             KfC(Ind3,Ind5,tM1)*UASf(All4,Ind4,bX)*
     &              UASfC(All3,3 + Ind5,tM1))*
     &           UASfC(All4,3 + Ind6,bX)))/(MW2*SW2)

	ENDLOOP(Ind1)
	ENDLOOP(Ind2)
	ENDLOOP(Ind3)
	ENDLOOP(Ind4)
	ENDLOOP(Ind5)
	ENDLOOP(Ind6)

	ENDLOOP(All3)
	ENDLOOP(All4)

#ifdef DETAILED_DEBUG
	DHIGGS "se_HmHp_nmfv =", se_HmHp_nmfv ENDL
#endif

	if( mssmpart .le. 2 ) return

	LOOP(Gen3, 1,3,1)

        se_HmHp_nmfv = se_HmHp_nmfv - 
     &    1/(4.D0*Pi)*(Alfa1L*TB2*
     &        (A0(Mf2(2,Gen3)) + k2*B1(k2,0.D0,Mf2(2,Gen3)))*
     &        Mf2(2,Gen3))/(MW2*SW2)

	ENDLOOP(Gen3)

	LOOP(All3, 1,6,1)

	tmp16 = A0(MASf2(All3,2))

	LOOP(Ind1, 1,3,1)

        se_HmHp_nmfv = se_HmHp_nmfv - 
     &    1/(16.D0*Pi)*(Alfa1L*tmp16*
     &        (C2B*MW2*UASf(All3,Ind1,2)*UASfC(All3,Ind1,2) - 
     &          2*(C2B*MW2*SW2 + CW2*TB2*Mf2(2,Ind1))*
     &           UASf(All3,3 + Ind1,2)*UASfC(All3,3 + Ind1,2)))/
     &      (CW2*MW2*SW2)

	ENDLOOP(Ind1)

	ENDLOOP(All3)

	LOOP(Gen3, 1,3,1)

	tmp17 = A0(MASf2(Gen3,1))

	LOOP(Ind1, 1,3,1)

        se_HmHp_nmfv = se_HmHp_nmfv - 
     &    1/(16.D0*Pi)*(Alfa1L*tmp17*
     &        (C2B*(MW2 - 2*CW2*MW2) - 2*CW2*TB2*Mf2(2,Ind1))*
     &        UASf(Gen3,Ind1,1)*UASfC(Gen3,Ind1,1))/(CW2*MW2*SW2)

	ENDLOOP(Ind1)

	ENDLOOP(Gen3)

	LOOP(Gen3, 1,3,1)
	LOOP(All3, 1,6,1)

	tmp18 = B0(k2,MASf2(All3,2),MASf2(Gen3,1))

	LOOP(Ind4, 1,3,1)
	LOOP(Ind3, 1,3,1)
	LOOP(Ind2, 1,3,1)
	LOOP(Ind1, 1,3,1)

        se_HmHp_nmfv = se_HmHp_nmfv + 
     &    1/(8.D0*Pi)*(Alfa1L*tmp18*
     &        (Delta(Ind1,Ind2)*(MW2*S2B - TB*Mf2(2,Ind1))*
     &           UASf(All3,Ind2,2) - 
     &          (TB*Kf(Ind1,Ind2,2) + 
     &             MUEC*Delta(Ind1,Ind2)*Mf(2,Ind1))*
     &           UASf(All3,3 + Ind2,2))*UASf(Gen3,Ind3,1)*
     &        (Delta(Ind3,Ind4)*(MW2*S2B - TB*Mf2(2,Ind3))*
     &           UASfC(All3,Ind4,2) - 
     &          (TB*KfC(Ind4,Ind3,2) + 
     &             MUE*Delta(Ind3,Ind4)*Mf(2,Ind3))*
     &           UASfC(All3,3 + Ind4,2))*UASfC(Gen3,Ind1,1))/
     &      (MW2*SW2)

	ENDLOOP(Ind1)
	ENDLOOP(Ind2)
	ENDLOOP(Ind3)
	ENDLOOP(Ind4)

	ENDLOOP(All3)
	ENDLOOP(Gen3)

#ifdef DETAILED_DEBUG
	DHIGGS "se_HmHp_nmfv =", se_HmHp_nmfv ENDL
#endif

	if( mssmpart .le. 3 ) return

        se_HmHp_nmfv = se_HmHp_nmfv - 
     &    1/(32.D0*Pi)*(Alfa1L*
     &        (2*((CBA2*CW2**2*(k2 + Mh02 - MW2) + 
     &                MW2*(2*CBA*CW2*S2B*SAB - S2B**2*SAB**2))*
     &              B0(k2,Mh02,MW2) - 
     &             MW2*((C2B*(C2B*SAB**2 + 4*CW2*SAB*SBA) + 
     &                   4*CW2**2*SBA2)*B0(k2,Mh02,MHp2) + 
     &                (4*CBA2*CW2**2 + 
     &                   C2B*(C2B*CAB**2 - 4*CAB*CBA*CW2))*
     &                 B0(k2,MHH2,MHp2)) - 
     &             (CAB*MW2*(CAB*S2B**2 - 2*CW2*S2B*SBA) - 
     &                CW2**2*(k2 + MHH2 - MW2)*SBA2)*
     &              B0(k2,MHH2,MW2)) - 
     &          CW2*((CW2*(1 - S2A*S2B) + C2A*C2B*SW2)*A0(Mh02) + 
     &             (CW2*(1 + S2A*S2B) - C2A*C2B*SW2)*A0(MHH2) + 
     &             C2B**2*(A0(MA02) + 4*A0(MHp2)) - 
     &             2*(C2B**2 - 6*CW2 - S2B**2)*A0(MW2) + 
     &             (CW2*(1 + 6*CW2 + S2B**2 - 12*SW2) - 
     &                C2B**2*SW2 + 6*SW2**2)*A0(MZ2) - 
     &             CW2*(2*(k2 + MA02 - MW2)*B0(k2,MA02,MW2) - 
     &                4*k2*
     &                 (B1(k2,MA02,MW2) + CBA2*B1(k2,Mh02,MW2) + 
     &                   SBA2*B1(k2,MHH2,MW2)) + 
     &                SW2*(8*MHp2*B0(k2,0.D0,MHp2) + 
     &                   k2*(8*B0(k2,0.D0,MHp2) - 16*B1(k2,MHp2,0.D0)))
     &                ) - (CW2 - SW2)**2*
     &              (2*MHp2*B0(k2,MHp2,MZ2) + 
     &                k2*(2*B0(k2,MHp2,MZ2) - 4*B1(k2,MHp2,MZ2)))))
     &        )/(CW2**2*SW2)

	LOOP(Cha3, 1,2,1)

	tmp19 = A0(MCha2(Cha3))

	LOOP(Neu3, 1,4,1)

	dup2 = (SW*ZNeu(Neu3,1) + CW*ZNeu(Neu3,2))/CW

	dup3 = (SW*ZNeuC(Neu3,1) + CW*ZNeuC(Neu3,2))/CW

        se_HmHp_nmfv = se_HmHp_nmfv - 
     &    1/(8.D0*Pi)*(Alfa1L*
     &        ((tmp19 - k2*B1(k2,MCha2(Cha3),MNeu2(Neu3)))*
     &           (4/sqrt2**2*
     &              (CB2*(sqrt2*(VCha(Cha3,1)*ZNeu(Neu3,4)) + 
     &                  dup2*VCha(Cha3,2))*
     &                (sqrt2*(VChaC(Cha3,1)*ZNeuC(Neu3,4)) + 
     &                  dup3*VChaC(Cha3,2))) + 
     &             SB2*(sqrt2*(dup2*UCha(Cha3,2)) - 
     &                2*UCha(Cha3,1)*ZNeu(Neu3,3))*
     &              (sqrt2*(dup3*UChaC(Cha3,2)) - 
     &                2*UChaC(Cha3,1)*ZNeuC(Neu3,3))) - 
     &          B0(k2,MCha2(Cha3),MNeu2(Neu3))*
     &           (4/sqrt2**2*
     &              (CB*SB*MCha(Cha3)*MNeu(Neu3)*
     &                ((-(sqrt2*(UCha(Cha3,1)*ZNeu(Neu3,3))) + 
     &                     dup2*UCha(Cha3,2))*
     &                   (sqrt2*(VCha(Cha3,1)*ZNeu(Neu3,4)) + 
     &                     dup2*VCha(Cha3,2)) + 
     &                  (-(sqrt2*(UChaC(Cha3,1)*ZNeuC(Neu3,3))) + 
     &                     dup3*UChaC(Cha3,2))*
     &                   (sqrt2*(VChaC(Cha3,1)*ZNeuC(Neu3,4)) + 
     &                     dup3*VChaC(Cha3,2)))) + 
     &             (k2 - MNeu2(Neu3))*
     &              (4/sqrt2**2*
     &                 (CB2*
     &                   (sqrt2*(VCha(Cha3,1)*ZNeu(Neu3,4)) + 
     &                     dup2*VCha(Cha3,2))*
     &                   (sqrt2*(VChaC(Cha3,1)*ZNeuC(Neu3,4)) + 
     &                     dup3*VChaC(Cha3,2))) + 
     &                SB2*(sqrt2*(dup2*UCha(Cha3,2)) - 
     &                   2*UCha(Cha3,1)*ZNeu(Neu3,3))*
     &                 (sqrt2*(dup3*UChaC(Cha3,2)) - 
     &                   2*UChaC(Cha3,1)*ZNeuC(Neu3,3))))))/SW2

	ENDLOOP(Neu3)

	ENDLOOP(Cha3)

#ifdef DETAILED_DEBUG
	DHIGGS "se_HmHp_nmfv =", se_HmHp_nmfv ENDL
#endif

	end

