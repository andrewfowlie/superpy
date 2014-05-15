* se_h0HH.F
* this file is part of FeynHiggs
* generated 05 Nov 2013 19:22

#include "externals.h"
#include "types.h"
#include "debug.h"


	ComplexType function se_h0HH(k2, bX)
	implicit none
	RealType k2
	integer bX


#include "FH.h"

	ComplexType se_h0HH_mfv, se_h0HH_nmfv
	external se_h0HH_mfv, se_h0HH_nmfv

	if( ibits(fv, 3, 2) .eq. 0 ) then
	  se_h0HH = se_h0HH_mfv(k2, bX)
	else
	  se_h0HH = se_h0HH_nmfv(k2, bX)
	endif

#ifdef DEBUG
	DHIGGS 'k2 =', k2 ENDL
	DHIGGS 'se_h0HH =', se_h0HH ENDL
#endif
	end


************************************************************************


	ComplexType function se_h0HH_mfv(k2,bX)
	implicit none
	integer bX
	RealType k2

#include "FH.h"
#include "looptools.h"

	integer Cha3, Cha4, Gen3, Neu3, Neu4, Sfe3, Sfe4

	ComplexType dup1, dup10, dup11, dup12, dup13, dup14, dup15
	ComplexType dup16, dup17, dup2, dup3, dup4, dup5, dup6, dup7
	ComplexType dup8, dup9, tmp1, tmp2

	se_h0HH_mfv = 0

	LOOP(Gen3, gM,3,1)

        se_h0HH_mfv = se_h0HH_mfv - 
     &    3/(4.D0*Pi)*(Alfa1L*CA*SA*Mf2(tM1,Gen3)*
     &        (A0(Mf2(tM1,Gen3)) + 
     &          k2*B1(k2,Mf2(tM1,Gen3),Mf2(tM1,Gen3)) + 
     &          2*B0(k2,Mf2(tM1,Gen3),Mf2(tM1,Gen3))*Mf2(tM1,Gen3))
     &        )/(MW2*SB2*SW2)

	ENDLOOP(Gen3)

	LOOP(Gen3, gM,3,1)
	LOOP(Sfe3, 1,2,1)

        se_h0HH_mfv = se_h0HH_mfv - 
     &    1/(16.D0*Pi)*(Alfa1L*S2A*A0(MSf2(Sfe3,tM1,Gen3))*
     &        ((MW2*SB2*(3 - 4*SW2) - 3*CW2*Mf2(tM1,Gen3))*
     &           USf2(Sfe3,1,tM1,Gen3) + 
     &          (4*MW2*SB2*SW2 - 3*CW2*Mf2(tM1,Gen3))*
     &           USf2(Sfe3,2,tM1,Gen3)))/(CW2*MW2*SB2*SW2)

	ENDLOOP(Sfe3)
	ENDLOOP(Gen3)

	LOOP(Gen3, gM,3,1)
	LOOP(Sfe4, 1,2,1)
	LOOP(Sfe3, 1,2,1)

        se_h0HH_mfv = se_h0HH_mfv - 
     &    1/(48.D0*Pi)*(Alfa1L*
     &        B0(k2,MSf2(Sfe3,tM1,Gen3),MSf2(Sfe4,tM1,Gen3))*
     &        (USf(Sfe4,1,tM1,Gen3)*
     &           ((MW*MZ*SAB*SB*(3 - 4*SW2) - 
     &                6*CA*CW*Mf2(tM1,Gen3))*USfC(Sfe3,1,tM1,Gen3)-
     &               3*CW*
     &              (CA*KfC(Gen3,Gen3,tM1) + MUE*SA*Mf(tM1,Gen3))*
     &              USfC(Sfe3,2,tM1,Gen3)) - 
     &          USf(Sfe4,2,tM1,Gen3)*
     &           (3*CW*(CA*Kf(Gen3,Gen3,tM1) + 
     &                MUEC*SA*Mf(tM1,Gen3))*USfC(Sfe3,1,tM1,Gen3)-
     &               2*(2*MW*MZ*SAB*SB*SW2 - 
     &                3*CA*CW*Mf2(tM1,Gen3))*USfC(Sfe3,2,tM1,Gen3))
     &          )*(USf(Sfe3,1,tM1,Gen3)*
     &           ((CAB*MW*MZ*SB*(3 - 4*SW2) + 
     &                6*CW*SA*Mf2(tM1,Gen3))*USfC(Sfe4,1,tM1,Gen3)+
     &               3*CW*
     &              (SA*KfC(Gen3,Gen3,tM1) - CA*MUE*Mf(tM1,Gen3))*
     &              USfC(Sfe4,2,tM1,Gen3)) + 
     &          USf(Sfe3,2,tM1,Gen3)*
     &           (3*CW*(SA*Kf(Gen3,Gen3,tM1) - 
     &                CA*MUEC*Mf(tM1,Gen3))*USfC(Sfe4,1,tM1,Gen3)+
     &               2*(2*CAB*MW*MZ*SB*SW2 + 
     &                3*CW*SA*Mf2(tM1,Gen3))*USfC(Sfe4,2,tM1,Gen3))
     &          ))/(CW2*MW2*SB2*SW2)

	ENDLOOP(Sfe3)
	ENDLOOP(Sfe4)
	ENDLOOP(Gen3)

#ifdef DETAILED_DEBUG
	DHIGGS "se_h0HH_mfv =", se_h0HH_mfv ENDL
#endif

	if( mssmpart .le. 1 ) return

	LOOP(Gen3, gM,3,1)

        se_h0HH_mfv = se_h0HH_mfv + 
     &    3/(4.D0*Pi)*(Alfa1L*CA*SA*Mf2(bX,Gen3)*
     &        (A0(Mf2(bX,Gen3)) + 
     &          k2*B1(k2,Mf2(bX,Gen3),Mf2(bX,Gen3)) + 
     &          2*B0(k2,Mf2(bX,Gen3),Mf2(bX,Gen3))*Mf2(bX,Gen3)))/
     &      (CB2*MW2*SW2)

	ENDLOOP(Gen3)

	LOOP(Gen3, gM,3,1)
	LOOP(Sfe3, 1,2,1)

        se_h0HH_mfv = se_h0HH_mfv + 
     &    1/(16.D0*Pi)*(Alfa1L*S2A*A0(MSf2(Sfe3,bX,Gen3))*
     &        ((CB2*MW2*(3 - 2*SW2) - 3*CW2*Mf2(bX,Gen3))*
     &           USf2(Sfe3,1,bX,Gen3) + 
     &          (2*CB2*MW2*SW2 - 3*CW2*Mf2(bX,Gen3))*
     &           USf2(Sfe3,2,bX,Gen3)))/(CB2*CW2*MW2*SW2)

	ENDLOOP(Sfe3)
	ENDLOOP(Gen3)

	LOOP(Gen3, gM,3,1)
	LOOP(Sfe4, 1,2,1)
	LOOP(Sfe3, 1,2,1)

        se_h0HH_mfv = se_h0HH_mfv - 
     &    1/(48.D0*Pi)*(Alfa1L*
     &        B0(k2,MSf2(Sfe3,bX,Gen3),MSf2(Sfe4,bX,Gen3))*
     &        (USf(Sfe4,1,bX,Gen3)*
     &           ((CB*MW*MZ*SAB*(3 - 2*SW2) - 
     &                6*CW*SA*Mf2(bX,Gen3))*USfC(Sfe3,1,bX,Gen3) - 
     &             3*CW*(SA*KfC(Gen3,Gen3,bX) + 
     &                CA*MUE*Mf(bX,Gen3))*USfC(Sfe3,2,bX,Gen3)) - 
     &          USf(Sfe4,2,bX,Gen3)*
     &           (3*CW*(SA*Kf(Gen3,Gen3,bX) + CA*MUEC*Mf(bX,Gen3))*
     &              USfC(Sfe3,1,bX,Gen3) - 
     &             2*(CB*MW*MZ*SAB*SW2 - 3*CW*SA*Mf2(bX,Gen3))*
     &              USfC(Sfe3,2,bX,Gen3)))*
     &        (USf(Sfe3,1,bX,Gen3)*
     &           ((CAB*CB*MW*MZ*(3 - 2*SW2) - 
     &                6*CA*CW*Mf2(bX,Gen3))*USfC(Sfe4,1,bX,Gen3) - 
     &             3*CW*(CA*KfC(Gen3,Gen3,bX) - 
     &                MUE*SA*Mf(bX,Gen3))*USfC(Sfe4,2,bX,Gen3)) - 
     &          USf(Sfe3,2,bX,Gen3)*
     &           (3*CW*(CA*Kf(Gen3,Gen3,bX) - MUEC*SA*Mf(bX,Gen3))*
     &              USfC(Sfe4,1,bX,Gen3) - 
     &             2*(CAB*CB*MW*MZ*SW2 - 3*CA*CW*Mf2(bX,Gen3))*
     &              USfC(Sfe4,2,bX,Gen3))))/(CB2*CW2*MW2*SW2)

	ENDLOOP(Sfe3)
	ENDLOOP(Sfe4)
	ENDLOOP(Gen3)

#ifdef DETAILED_DEBUG
	DHIGGS "se_h0HH_mfv =", se_h0HH_mfv ENDL
#endif

	if( mssmpart .le. 2 ) return

	LOOP(Gen3, 1,3,1)

        se_h0HH_mfv = se_h0HH_mfv - 
     &    1/(16.D0*Pi)*(Alfa1L*
     &        (CB2*MW2*(S2A*A0(MSf2(1,1,Gen3)) + 
     &             CAB*MZ2*SAB*B0(k2,MSf2(1,1,Gen3),MSf2(1,1,Gen3))
     &             ) - 4*CA*CW2*SA*Mf2(2,Gen3)*
     &           (A0(Mf2(2,Gen3)) + 
     &             k2*B1(k2,Mf2(2,Gen3),Mf2(2,Gen3)) + 
     &             2*B0(k2,Mf2(2,Gen3),Mf2(2,Gen3))*Mf2(2,Gen3))))/
     &      (CB2*CW2*MW2*SW2)

	ENDLOOP(Gen3)

	LOOP(Sfe3, 1,2,1)
	LOOP(Gen3, 1,3,1)

        se_h0HH_mfv = se_h0HH_mfv + 
     &    1/(16.D0*Pi)*(Alfa1L*S2A*A0(MSf2(Sfe3,2,Gen3))*
     &        ((CB2*MW2*(1 - 2*SW2) - CW2*Mf2(2,Gen3))*
     &           USf2(Sfe3,1,2,Gen3) + 
     &          (2*CB2*MW2*SW2 - CW2*Mf2(2,Gen3))*
     &           USf2(Sfe3,2,2,Gen3)))/(CB2*CW2*MW2*SW2)

	ENDLOOP(Gen3)
	ENDLOOP(Sfe3)

	LOOP(Sfe4, 1,2,1)
	LOOP(Sfe3, 1,2,1)
	LOOP(Gen3, 1,3,1)

        se_h0HH_mfv = se_h0HH_mfv - 
     &    1/(16.D0*Pi)*(Alfa1L*
     &        B0(k2,MSf2(Sfe3,2,Gen3),MSf2(Sfe4,2,Gen3))*
     &        (USf(Sfe4,1,2,Gen3)*
     &           ((CB*MW*MZ*SAB*(1 - 2*SW2) - 2*CW*SA*Mf2(2,Gen3))*
     &              USfC(Sfe3,1,2,Gen3) - 
     &             CW*(SA*KfC(Gen3,Gen3,2) + CA*MUE*Mf(2,Gen3))*
     &              USfC(Sfe3,2,2,Gen3)) - 
     &          USf(Sfe4,2,2,Gen3)*
     &           (CW*(SA*Kf(Gen3,Gen3,2) + CA*MUEC*Mf(2,Gen3))*
     &              USfC(Sfe3,1,2,Gen3) - 
     &             2*(CB*MW*MZ*SAB*SW2 - CW*SA*Mf2(2,Gen3))*
     &              USfC(Sfe3,2,2,Gen3)))*
     &        (USf(Sfe3,1,2,Gen3)*
     &           ((CAB*CB*MW*MZ*(1 - 2*SW2) - 2*CA*CW*Mf2(2,Gen3))*
     &              USfC(Sfe4,1,2,Gen3) - 
     &             CW*(CA*KfC(Gen3,Gen3,2) - MUE*SA*Mf(2,Gen3))*
     &              USfC(Sfe4,2,2,Gen3)) - 
     &          USf(Sfe3,2,2,Gen3)*
     &           (CW*(CA*Kf(Gen3,Gen3,2) - MUEC*SA*Mf(2,Gen3))*
     &              USfC(Sfe4,1,2,Gen3) - 
     &             2*(CAB*CB*MW*MZ*SW2 - CA*CW*Mf2(2,Gen3))*
     &              USfC(Sfe4,2,2,Gen3))))/(CB2*CW2*MW2*SW2)

	ENDLOOP(Gen3)
	ENDLOOP(Sfe3)
	ENDLOOP(Sfe4)

#ifdef DETAILED_DEBUG
	DHIGGS "se_h0HH_mfv =", se_h0HH_mfv ENDL
#endif

	if( mssmpart .le. 3 ) return

        se_h0HH_mfv = se_h0HH_mfv - 
     &    1/(32.D0*Pi)*(Alfa1L*
     &        (MW2*(SAB*(-6*C2A*S2A*SAB*B0(k2,Mh02,Mh02) + 
     &                CAB*(C2B**2*B0(k2,MA02,MA02) + 
     &                   3*C2A**2*B0(k2,Mh02,Mh02))) + 
     &             3*C2A*CAB*(2*CAB*S2A + C2A*SAB)*B0(k2,MHH2,MHH2)
     &             ) + 2*((CAB*MW2*S2B**2*SAB - 
     &                CBA*CW2*(k2 + MA02)*SBA)*B0(k2,MA02,MZ2) + 
     &             MW2*((4*CAB*S2A**2*SAB - 
     &                   C2A*
     &                    (2*CAB**2*S2A + C2A*CAB*SAB - 
     &                      2*S2A*SAB**2))*B0(k2,Mh02,MHH2) + 
     &                (-4*CBA*CW2**2*SBA + 
     &                   C2B*
     &                    ((C2B*CAB - 2*CBA*CW2)*SAB + 
     &                      2*CAB*CW2*SBA))*B0(k2,MHp2,MHp2)) + 
     &             (C2B**2*CAB*MW2*SAB + 
     &                CBA*CW2**2*(2*k2 - 12*MW2)*SBA)*
     &              B0(k2,MW2,MW2)) + 
     &          (C2B**2*CAB*MW2*SAB - 
     &             2*CBA*(7*MW2 - CW2*(k2 + MZ2))*SBA)*
     &           B0(k2,MZ2,MZ2) - 
     &          CW2*((C2A*CW2*S2B + C2B*S2A*SW2)*
     &              (2*A0(MHp2) - 2*A0(MW2)) + 
     &             S2A*(C2A*(3*A0(Mh02) - 3*A0(MHH2)) + 
     &                C2B*(A0(MA02) - A0(MZ2))) - 
     &             4*CBA*k2*SBA*B1(k2,MA02,MZ2)) + 
     &          CBA*CW2**2*k2*SBA*
     &           (8*B1(k2,MHp2,MW2) - 8*B1(k2,MW2,MW2)) - 
     &          4*((CAB*MW2*S2B*(CBA*CW2 - S2B*SAB) + 
     &                CW2*(CBA*CW2*(k2 + MHp2 - MW2) + 
     &                   MW2*S2B*SAB)*SBA)*B0(k2,MHp2,MW2) + 
     &             CBA*CW2*k2*SBA*B1(k2,MZ2,MZ2))))/(CW2**2*SW2)

	LOOP(Cha4, 1,2,1)

	tmp1 = A0(MCha2(Cha4))

	LOOP(Cha3, 1,2,1)

        dup1 = SA*UCha(Cha4,2)*VCha(Cha3,1) - 
     &    CA*UCha(Cha4,1)*VCha(Cha3,2)

        dup2 = CA*UCha(Cha3,2)*VCha(Cha4,1) + 
     &    SA*UCha(Cha3,1)*VCha(Cha4,2)

        dup3 = CA*UChaC(Cha4,2)*VChaC(Cha3,1) + 
     &    SA*UChaC(Cha4,1)*VChaC(Cha3,2)

        dup4 = SA*UChaC(Cha3,2)*VChaC(Cha4,1) - 
     &    CA*UChaC(Cha3,1)*VChaC(Cha4,2)

        se_h0HH_mfv = se_h0HH_mfv + 
     &    1/(4.D0*Pi)*(Alfa1L*
     &        ((dup1*dup3 + dup2*dup4)*
     &           (tmp1 + k2*B1(k2,MCha2(Cha3),MCha2(Cha4))) + 
     &          B0(k2,MCha2(Cha3),MCha2(Cha4))*
     &           (dup4*(dup3*MCha(Cha3)*MCha(Cha4) + 
     &                dup2*MCha2(Cha3)) + 
     &             dup1*(dup2*MCha(Cha3)*MCha(Cha4) + 
     &                dup3*MCha2(Cha3)))))/SW2

	ENDLOOP(Cha3)

	ENDLOOP(Cha4)

	LOOP(Neu4, 1,4,1)

	tmp2 = A0(MNeu2(Neu4))

	LOOP(Neu3, 1,4,1)

	dup5 = SW*ZNeu(Neu3,1) - CW*ZNeu(Neu3,2)

	dup6 = SA*ZNeu(Neu3,3) + CA*ZNeu(Neu3,4)

	dup7 = CA*ZNeu(Neu3,3) - SA*ZNeu(Neu3,4)

	dup8 = SW*ZNeu(Neu4,1) - CW*ZNeu(Neu4,2)

	dup9 = SA*ZNeu(Neu4,3) + CA*ZNeu(Neu4,4)

	dup10 = CA*ZNeu(Neu4,3) - SA*ZNeu(Neu4,4)

	dup11 = SW*ZNeuC(Neu3,1) - CW*ZNeuC(Neu3,2)

	dup12 = SA*ZNeuC(Neu3,3) + CA*ZNeuC(Neu3,4)

	dup13 = CA*ZNeuC(Neu3,3) - SA*ZNeuC(Neu3,4)

	dup14 = SW*ZNeuC(Neu4,1) - CW*ZNeuC(Neu4,2)

	dup15 = SA*ZNeuC(Neu4,3) + CA*ZNeuC(Neu4,4)

	dup16 = CA*ZNeuC(Neu4,3) - SA*ZNeuC(Neu4,4)

        dup17 = dup10*(dup12*dup14 + dup11*dup15)*dup5 + 
     &    dup12*dup14*dup7*dup8 + 
     &    dup13*dup14*(dup6*dup8 + dup5*dup9) + 
     &    dup11*((dup16*dup6 + dup15*dup7)*dup8 + dup16*dup5*dup9)

        se_h0HH_mfv = se_h0HH_mfv + 
     &    1/(16.D0*Pi)*(Alfa1L*
     &        (dup17*(tmp2 + k2*B1(k2,MNeu2(Neu3),MNeu2(Neu4))) + 
     &          B0(k2,MNeu2(Neu3),MNeu2(Neu4))*
     &           ((dup15*(dup11*dup13*dup14 + dup11**2*dup16) + 
     &                dup12*(dup13*dup14**2 + dup11*dup14*dup16) + 
     &                dup10*dup5*(dup6*dup8 + dup5*dup9) + 
     &                dup7*(dup6*dup8**2 + dup5*dup8*dup9))*
     &              MNeu(Neu3)*MNeu(Neu4) + dup17*MNeu2(Neu3))))/
     &      (CW2*SW2)

	ENDLOOP(Neu3)

	ENDLOOP(Neu4)

#ifdef DETAILED_DEBUG
	DHIGGS "se_h0HH_mfv =", se_h0HH_mfv ENDL
#endif

	end


************************************************************************


	ComplexType function se_h0HH_nmfv(k2,bX)
	implicit none
	integer bX
	RealType k2

#include "FH.h"
#include "looptools.h"

	integer All3, All4, Cha3, Cha4, Gen3, Ind1, Ind2, Ind3, Ind4
	integer Neu3, Neu4

	ComplexType dup1, dup10, dup11, dup12, dup13, dup14, dup15
	ComplexType dup16, dup17, dup2, dup3, dup4, dup5, dup6, dup7
	ComplexType dup8, dup9, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6
	ComplexType tmp7, tmp8

	se_h0HH_nmfv = 0

	LOOP(Gen3, gM,3,1)

        se_h0HH_nmfv = se_h0HH_nmfv - 
     &    3/(4.D0*Pi)*(Alfa1L*CA*SA*Mf2(tM1,Gen3)*
     &        (A0(Mf2(tM1,Gen3)) + 
     &          k2*B1(k2,Mf2(tM1,Gen3),Mf2(tM1,Gen3)) + 
     &          2*B0(k2,Mf2(tM1,Gen3),Mf2(tM1,Gen3))*Mf2(tM1,Gen3))
     &        )/(MW2*SB2*SW2)

	ENDLOOP(Gen3)

	LOOP(All3, 1,6,gM)

	tmp1 = A0(MASf2(All3,tM1))

	LOOP(Ind1, 1,3,1)

        se_h0HH_nmfv = se_h0HH_nmfv - 
     &    1/(16.D0*Pi)*(Alfa1L*S2A*tmp1*
     &        ((MW2*SB2*(3 - 4*SW2) - 3*CW2*Mf2(tM1,Ind1))*
     &           UASf(All3,Ind1,tM1)*UASfC(All3,Ind1,tM1) + 
     &          (4*MW2*SB2*SW2 - 3*CW2*Mf2(tM1,Ind1))*
     &           UASf(All3,3 + Ind1,tM1)*UASfC(All3,3 + Ind1,tM1)))
     &       /(CW2*MW2*SB2*SW2)

	ENDLOOP(Ind1)

	ENDLOOP(All3)

	LOOP(All4, 1,6,gM)
	LOOP(All3, 1,6,gM)

	tmp2 = B0(k2,MASf2(All3,tM1),MASf2(All4,tM1))

	LOOP(Ind4, 1,3,1)
	LOOP(Ind3, 1,3,1)
	LOOP(Ind2, 1,3,1)
	LOOP(Ind1, 1,3,1)

        se_h0HH_nmfv = se_h0HH_nmfv - 
     &    1/(48.D0*Pi)*(Alfa1L*tmp2*
     &        ((Delta(Ind1,Ind2)*
     &              (MW*MZ*SAB*SB*(3 - 4*SW2) - 
     &                6*CA*CW*Mf2(tM1,Ind1))*UASf(All4,Ind2,tM1) - 
     &             3*CW*(CA*Kf(Ind1,Ind2,tM1) + 
     &                MUEC*SA*Delta(Ind1,Ind2)*Mf(tM1,Ind1))*
     &              UASf(All4,3 + Ind2,tM1))*UASfC(All3,Ind1,tM1)-
     &            (3*CW*(CA*KfC(Ind2,Ind1,tM1) + 
     &                MUE*SA*Delta(Ind1,Ind2)*Mf(tM1,Ind1))*
     &              UASf(All4,Ind2,tM1) - 
     &             2*Delta(Ind1,Ind2)*
     &              (2*MW*MZ*SAB*SB*SW2 - 3*CA*CW*Mf2(tM1,Ind1))*
     &              UASf(All4,3 + Ind2,tM1))*
     &           UASfC(All3,3 + Ind1,tM1))*
     &        ((Delta(Ind3,Ind4)*
     &              (CAB*MW*MZ*SB*(3 - 4*SW2) + 
     &                6*CW*SA*Mf2(tM1,Ind3))*UASf(All3,Ind4,tM1) + 
     &             3*CW*(SA*Kf(Ind3,Ind4,tM1) - 
     &                CA*MUEC*Delta(Ind3,Ind4)*Mf(tM1,Ind3))*
     &              UASf(All3,3 + Ind4,tM1))*UASfC(All4,Ind3,tM1)+
     &            (3*CW*(SA*KfC(Ind4,Ind3,tM1) - 
     &                CA*MUE*Delta(Ind3,Ind4)*Mf(tM1,Ind3))*
     &              UASf(All3,Ind4,tM1) + 
     &             2*Delta(Ind3,Ind4)*
     &              (2*CAB*MW*MZ*SB*SW2 + 3*CW*SA*Mf2(tM1,Ind3))*
     &              UASf(All3,3 + Ind4,tM1))*
     &           UASfC(All4,3 + Ind3,tM1)))/(CW2*MW2*SB2*SW2)

	ENDLOOP(Ind1)
	ENDLOOP(Ind2)
	ENDLOOP(Ind3)
	ENDLOOP(Ind4)

	ENDLOOP(All3)
	ENDLOOP(All4)

#ifdef DETAILED_DEBUG
	DHIGGS "se_h0HH_nmfv =", se_h0HH_nmfv ENDL
#endif

	if( mssmpart .le. 1 ) return

	LOOP(Gen3, gM,3,1)

        se_h0HH_nmfv = se_h0HH_nmfv + 
     &    3/(4.D0*Pi)*(Alfa1L*CA*SA*Mf2(bX,Gen3)*
     &        (A0(Mf2(bX,Gen3)) + 
     &          k2*B1(k2,Mf2(bX,Gen3),Mf2(bX,Gen3)) + 
     &          2*B0(k2,Mf2(bX,Gen3),Mf2(bX,Gen3))*Mf2(bX,Gen3)))/
     &      (CB2*MW2*SW2)

	ENDLOOP(Gen3)

	LOOP(All3, 1,6,gM)

	tmp3 = A0(MASf2(All3,bX))

	LOOP(Ind1, 1,3,1)

        se_h0HH_nmfv = se_h0HH_nmfv + 
     &    1/(16.D0*Pi)*(Alfa1L*S2A*tmp3*
     &        ((CB2*MW2*(3 - 2*SW2) - 3*CW2*Mf2(bX,Ind1))*
     &           UASf(All3,Ind1,bX)*UASfC(All3,Ind1,bX) + 
     &          (2*CB2*MW2*SW2 - 3*CW2*Mf2(bX,Ind1))*
     &           UASf(All3,3 + Ind1,bX)*UASfC(All3,3 + Ind1,bX)))/
     &      (CB2*CW2*MW2*SW2)

	ENDLOOP(Ind1)

	ENDLOOP(All3)

	LOOP(All4, 1,6,gM)
	LOOP(All3, 1,6,gM)

	tmp4 = B0(k2,MASf2(All3,bX),MASf2(All4,bX))

	LOOP(Ind4, 1,3,1)
	LOOP(Ind3, 1,3,1)
	LOOP(Ind2, 1,3,1)
	LOOP(Ind1, 1,3,1)

        se_h0HH_nmfv = se_h0HH_nmfv - 
     &    1/(48.D0*Pi)*(Alfa1L*tmp4*
     &        ((Delta(Ind3,Ind4)*
     &              (CB*MW*MZ*SAB*(3 - 2*SW2) - 
     &                6*CW*SA*Mf2(bX,Ind3))*UASf(All4,Ind4,bX) - 
     &             3*CW*(SA*Kf(Ind3,Ind4,bX) + 
     &                CA*MUEC*Delta(Ind3,Ind4)*Mf(bX,Ind3))*
     &              UASf(All4,3 + Ind4,bX))*UASfC(All3,Ind3,bX) - 
     &          (3*CW*(SA*KfC(Ind4,Ind3,bX) + 
     &                CA*MUE*Delta(Ind3,Ind4)*Mf(bX,Ind3))*
     &              UASf(All4,Ind4,bX) - 
     &             2*Delta(Ind3,Ind4)*
     &              (CB*MW*MZ*SAB*SW2 - 3*CW*SA*Mf2(bX,Ind3))*
     &              UASf(All4,3 + Ind4,bX))*UASfC(All3,3 + Ind3,bX)
     &          )*((Delta(Ind1,Ind2)*
     &              (CAB*CB*MW*MZ*(3 - 2*SW2) - 
     &                6*CA*CW*Mf2(bX,Ind1))*UASf(All3,Ind2,bX) - 
     &             3*CW*(CA*Kf(Ind1,Ind2,bX) - 
     &                MUEC*SA*Delta(Ind1,Ind2)*Mf(bX,Ind1))*
     &              UASf(All3,3 + Ind2,bX))*UASfC(All4,Ind1,bX) - 
     &          (3*CW*(CA*KfC(Ind2,Ind1,bX) - 
     &                MUE*SA*Delta(Ind1,Ind2)*Mf(bX,Ind1))*
     &              UASf(All3,Ind2,bX) - 
     &             2*Delta(Ind1,Ind2)*
     &              (CAB*CB*MW*MZ*SW2 - 3*CA*CW*Mf2(bX,Ind1))*
     &              UASf(All3,3 + Ind2,bX))*UASfC(All4,3 + Ind1,bX)
     &          ))/(CB2*CW2*MW2*SW2)

	ENDLOOP(Ind1)
	ENDLOOP(Ind2)
	ENDLOOP(Ind3)
	ENDLOOP(Ind4)

	ENDLOOP(All3)
	ENDLOOP(All4)

#ifdef DETAILED_DEBUG
	DHIGGS "se_h0HH_nmfv =", se_h0HH_nmfv ENDL
#endif

	if( mssmpart .le. 2 ) return

	LOOP(Gen3, 1,3,1)

        se_h0HH_nmfv = se_h0HH_nmfv - 
     &    1/(16.D0*Pi)*(Alfa1L*
     &        (CB2*MW2*(S2A*A0(MASf2(Gen3,1)) + 
     &             CAB*MZ2*SAB*B0(k2,MASf2(Gen3,1),MASf2(Gen3,1)))-
     &            4*CA*CW2*SA*Mf2(2,Gen3)*
     &           (A0(Mf2(2,Gen3)) + 
     &             k2*B1(k2,Mf2(2,Gen3),Mf2(2,Gen3)) + 
     &             2*B0(k2,Mf2(2,Gen3),Mf2(2,Gen3))*Mf2(2,Gen3))))/
     &      (CB2*CW2*MW2*SW2)

	ENDLOOP(Gen3)

	LOOP(All3, 1,6,1)

	tmp5 = A0(MASf2(All3,2))

	LOOP(Ind1, 1,3,1)

        se_h0HH_nmfv = se_h0HH_nmfv + 
     &    1/(16.D0*Pi)*(Alfa1L*S2A*tmp5*
     &        ((CB2*MW2*(1 - 2*SW2) - CW2*Mf2(2,Ind1))*
     &           UASf(All3,Ind1,2)*UASfC(All3,Ind1,2) + 
     &          (2*CB2*MW2*SW2 - CW2*Mf2(2,Ind1))*
     &           UASf(All3,3 + Ind1,2)*UASfC(All3,3 + Ind1,2)))/
     &      (CB2*CW2*MW2*SW2)

	ENDLOOP(Ind1)

	ENDLOOP(All3)

	LOOP(All4, 1,6,1)
	LOOP(All3, 1,6,1)

	tmp6 = B0(k2,MASf2(All3,2),MASf2(All4,2))

	LOOP(Ind4, 1,3,1)
	LOOP(Ind3, 1,3,1)
	LOOP(Ind2, 1,3,1)
	LOOP(Ind1, 1,3,1)

        se_h0HH_nmfv = se_h0HH_nmfv - 
     &    1/(16.D0*Pi)*(Alfa1L*tmp6*
     &        ((Delta(Ind1,Ind2)*
     &              (CB*MW*MZ*SAB*(1 - 2*SW2) - 
     &                2*CW*SA*Mf2(2,Ind1))*UASf(All4,Ind2,2) - 
     &             CW*(SA*Kf(Ind1,Ind2,2) + 
     &                CA*MUEC*Delta(Ind1,Ind2)*Mf(2,Ind1))*
     &              UASf(All4,3 + Ind2,2))*UASfC(All3,Ind1,2) - 
     &          (CW*(SA*KfC(Ind2,Ind1,2) + 
     &                CA*MUE*Delta(Ind1,Ind2)*Mf(2,Ind1))*
     &              UASf(All4,Ind2,2) - 
     &             2*Delta(Ind1,Ind2)*
     &              (CB*MW*MZ*SAB*SW2 - CW*SA*Mf2(2,Ind1))*
     &              UASf(All4,3 + Ind2,2))*UASfC(All3,3 + Ind1,2))*
     &        ((Delta(Ind3,Ind4)*
     &              (CAB*CB*MW*MZ*(1 - 2*SW2) - 
     &                2*CA*CW*Mf2(2,Ind3))*UASf(All3,Ind4,2) - 
     &             CW*(CA*Kf(Ind3,Ind4,2) - 
     &                MUEC*SA*Delta(Ind3,Ind4)*Mf(2,Ind3))*
     &              UASf(All3,3 + Ind4,2))*UASfC(All4,Ind3,2) - 
     &          (CW*(CA*KfC(Ind4,Ind3,2) - 
     &                MUE*SA*Delta(Ind3,Ind4)*Mf(2,Ind3))*
     &              UASf(All3,Ind4,2) - 
     &             2*Delta(Ind3,Ind4)*
     &              (CAB*CB*MW*MZ*SW2 - CA*CW*Mf2(2,Ind3))*
     &              UASf(All3,3 + Ind4,2))*UASfC(All4,3 + Ind3,2)))
     &       /(CB2*CW2*MW2*SW2)

	ENDLOOP(Ind1)
	ENDLOOP(Ind2)
	ENDLOOP(Ind3)
	ENDLOOP(Ind4)

	ENDLOOP(All3)
	ENDLOOP(All4)

#ifdef DETAILED_DEBUG
	DHIGGS "se_h0HH_nmfv =", se_h0HH_nmfv ENDL
#endif

	if( mssmpart .le. 3 ) return

        se_h0HH_nmfv = se_h0HH_nmfv - 
     &    1/(32.D0*Pi)*(Alfa1L*
     &        (MW2*(SAB*(-6*C2A*S2A*SAB*B0(k2,Mh02,Mh02) + 
     &                CAB*(C2B**2*B0(k2,MA02,MA02) + 
     &                   3*C2A**2*B0(k2,Mh02,Mh02))) + 
     &             3*C2A*CAB*(2*CAB*S2A + C2A*SAB)*B0(k2,MHH2,MHH2)
     &             ) + 2*((CAB*MW2*S2B**2*SAB - 
     &                CBA*CW2*(k2 + MA02)*SBA)*B0(k2,MA02,MZ2) + 
     &             MW2*((4*CAB*S2A**2*SAB - 
     &                   C2A*
     &                    (2*CAB**2*S2A + C2A*CAB*SAB - 
     &                      2*S2A*SAB**2))*B0(k2,Mh02,MHH2) + 
     &                (-4*CBA*CW2**2*SBA + 
     &                   C2B*
     &                    ((C2B*CAB - 2*CBA*CW2)*SAB + 
     &                      2*CAB*CW2*SBA))*B0(k2,MHp2,MHp2)) + 
     &             (C2B**2*CAB*MW2*SAB + 
     &                CBA*CW2**2*(2*k2 - 12*MW2)*SBA)*
     &              B0(k2,MW2,MW2)) + 
     &          (C2B**2*CAB*MW2*SAB - 
     &             2*CBA*(7*MW2 - CW2*(k2 + MZ2))*SBA)*
     &           B0(k2,MZ2,MZ2) - 
     &          CW2*((C2A*CW2*S2B + C2B*S2A*SW2)*
     &              (2*A0(MHp2) - 2*A0(MW2)) + 
     &             S2A*(C2A*(3*A0(Mh02) - 3*A0(MHH2)) + 
     &                C2B*(A0(MA02) - A0(MZ2))) - 
     &             4*CBA*k2*SBA*B1(k2,MA02,MZ2)) + 
     &          CBA*CW2**2*k2*SBA*
     &           (8*B1(k2,MHp2,MW2) - 8*B1(k2,MW2,MW2)) - 
     &          4*((CAB*MW2*S2B*(CBA*CW2 - S2B*SAB) + 
     &                CW2*(CBA*CW2*(k2 + MHp2 - MW2) + 
     &                   MW2*S2B*SAB)*SBA)*B0(k2,MHp2,MW2) + 
     &             CBA*CW2*k2*SBA*B1(k2,MZ2,MZ2))))/(CW2**2*SW2)

	LOOP(Cha4, 1,2,1)

	tmp7 = A0(MCha2(Cha4))

	LOOP(Cha3, 1,2,1)

        dup1 = SA*UCha(Cha4,2)*VCha(Cha3,1) - 
     &    CA*UCha(Cha4,1)*VCha(Cha3,2)

        dup2 = CA*UCha(Cha3,2)*VCha(Cha4,1) + 
     &    SA*UCha(Cha3,1)*VCha(Cha4,2)

        dup3 = CA*UChaC(Cha4,2)*VChaC(Cha3,1) + 
     &    SA*UChaC(Cha4,1)*VChaC(Cha3,2)

        dup4 = SA*UChaC(Cha3,2)*VChaC(Cha4,1) - 
     &    CA*UChaC(Cha3,1)*VChaC(Cha4,2)

        se_h0HH_nmfv = se_h0HH_nmfv + 
     &    1/(4.D0*Pi)*(Alfa1L*
     &        ((dup1*dup3 + dup2*dup4)*
     &           (tmp7 + k2*B1(k2,MCha2(Cha3),MCha2(Cha4))) + 
     &          B0(k2,MCha2(Cha3),MCha2(Cha4))*
     &           (dup4*(dup3*MCha(Cha3)*MCha(Cha4) + 
     &                dup2*MCha2(Cha3)) + 
     &             dup1*(dup2*MCha(Cha3)*MCha(Cha4) + 
     &                dup3*MCha2(Cha3)))))/SW2

	ENDLOOP(Cha3)

	ENDLOOP(Cha4)

	LOOP(Neu4, 1,4,1)

	tmp8 = A0(MNeu2(Neu4))

	LOOP(Neu3, 1,4,1)

	dup5 = SW*ZNeu(Neu3,1) - CW*ZNeu(Neu3,2)

	dup6 = SA*ZNeu(Neu3,3) + CA*ZNeu(Neu3,4)

	dup7 = CA*ZNeu(Neu3,3) - SA*ZNeu(Neu3,4)

	dup8 = SW*ZNeu(Neu4,1) - CW*ZNeu(Neu4,2)

	dup9 = SA*ZNeu(Neu4,3) + CA*ZNeu(Neu4,4)

	dup10 = CA*ZNeu(Neu4,3) - SA*ZNeu(Neu4,4)

	dup11 = SW*ZNeuC(Neu3,1) - CW*ZNeuC(Neu3,2)

	dup12 = SA*ZNeuC(Neu3,3) + CA*ZNeuC(Neu3,4)

	dup13 = CA*ZNeuC(Neu3,3) - SA*ZNeuC(Neu3,4)

	dup14 = SW*ZNeuC(Neu4,1) - CW*ZNeuC(Neu4,2)

	dup15 = SA*ZNeuC(Neu4,3) + CA*ZNeuC(Neu4,4)

	dup16 = CA*ZNeuC(Neu4,3) - SA*ZNeuC(Neu4,4)

        dup17 = dup10*(dup12*dup14 + dup11*dup15)*dup5 + 
     &    dup12*dup14*dup7*dup8 + 
     &    dup13*dup14*(dup6*dup8 + dup5*dup9) + 
     &    dup11*((dup16*dup6 + dup15*dup7)*dup8 + dup16*dup5*dup9)

        se_h0HH_nmfv = se_h0HH_nmfv + 
     &    1/(16.D0*Pi)*(Alfa1L*
     &        (dup17*(tmp8 + k2*B1(k2,MNeu2(Neu3),MNeu2(Neu4))) + 
     &          B0(k2,MNeu2(Neu3),MNeu2(Neu4))*
     &           ((dup15*(dup11*dup13*dup14 + dup11**2*dup16) + 
     &                dup12*(dup13*dup14**2 + dup11*dup14*dup16) + 
     &                dup10*dup5*(dup6*dup8 + dup5*dup9) + 
     &                dup7*(dup6*dup8**2 + dup5*dup8*dup9))*
     &              MNeu(Neu3)*MNeu(Neu4) + dup17*MNeu2(Neu3))))/
     &      (CW2*SW2)

	ENDLOOP(Neu3)

	ENDLOOP(Neu4)

#ifdef DETAILED_DEBUG
	DHIGGS "se_h0HH_nmfv =", se_h0HH_nmfv ENDL
#endif

	end

