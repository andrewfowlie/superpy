* dse_HHHH.F
* this file is part of FeynHiggs
* generated 05 Nov 2013 19:22

#include "externals.h"
#include "types.h"
#include "debug.h"


	ComplexType function dse_HHHH(k2, bX)
	implicit none
	RealType k2
	integer bX


#include "FH.h"

	ComplexType dse_HHHH_mfv, dse_HHHH_nmfv
	external dse_HHHH_mfv, dse_HHHH_nmfv

	if( ibits(fv, 3, 2) .eq. 0 ) then
	  dse_HHHH = dse_HHHH_mfv(k2, bX)
	else
	  dse_HHHH = dse_HHHH_nmfv(k2, bX)
	endif

#ifdef DEBUG
	DHIGGS 'k2 =', k2 ENDL
	DHIGGS 'dse_HHHH =', dse_HHHH ENDL
#endif
	end


************************************************************************


	ComplexType function dse_HHHH_mfv(k2,bX)
	implicit none
	integer bX
	RealType k2

#include "FH.h"
#include "looptools.h"

	integer Cha3, Cha4, Gen3, Neu3, Neu4, Sfe3, Sfe4

	ComplexType dup1, dup10, dup11, dup12, dup13, dup14, dup15
	ComplexType dup16, dup17, dup18, dup19, dup2, dup20, dup21
	ComplexType dup22, dup23, dup24, dup3, dup4, dup5, dup6
	ComplexType dup7, dup8, dup9

	dse_HHHH_mfv = 0

	LOOP(Gen3, gM,3,1)

        dse_HHHH_mfv = dse_HHHH_mfv - 
     &    3/(4.D0*Pi)*(Alfa1L*SA2*Mf2(tM1,Gen3)*
     &        (B1(k2,Mf2(tM1,Gen3),Mf2(tM1,Gen3)) + 
     &          k2*DB1(k2,Mf2(tM1,Gen3),Mf2(tM1,Gen3)) + 
     &          2*DB0(k2,Mf2(tM1,Gen3),Mf2(tM1,Gen3))*Mf2(tM1,Gen3)
     &          ))/(MW2*SB2*SW2)

	ENDLOOP(Gen3)

	LOOP(Gen3, gM,3,1)
	LOOP(Sfe4, 1,2,1)
	LOOP(Sfe3, 1,2,1)

	dup1 = SA*KfC(Gen3,Gen3,tM1) - CA*MUE*Mf(tM1,Gen3)

	dup2 = SA*Kf(Gen3,Gen3,tM1) - CA*MUEC*Mf(tM1,Gen3)

	dup3 = 2*CAB*MW*MZ*SB*SW2 + 3*CW*SA*Mf2(tM1,Gen3)

	dup4 = CAB*MW*MZ*SB*(3 - 4*SW2) + 6*CW*SA*Mf2(tM1,Gen3)

        dse_HHHH_mfv = dse_HHHH_mfv + 
     &    1/(48.D0*Pi)*(Alfa1L*
     &        DB0(k2,MSf2(Sfe3,tM1,Gen3),MSf2(Sfe4,tM1,Gen3))*
     &        (USf(Sfe4,1,tM1,Gen3)*
     &           (dup4*USfC(Sfe3,1,tM1,Gen3) + 
     &             3*CW*dup1*USfC(Sfe3,2,tM1,Gen3)) + 
     &          USf(Sfe4,2,tM1,Gen3)*
     &           (3*CW*dup2*USfC(Sfe3,1,tM1,Gen3) + 
     &             2*dup3*USfC(Sfe3,2,tM1,Gen3)))*
     &        (USf(Sfe3,1,tM1,Gen3)*
     &           (dup4*USfC(Sfe4,1,tM1,Gen3) + 
     &             3*CW*dup1*USfC(Sfe4,2,tM1,Gen3)) + 
     &          USf(Sfe3,2,tM1,Gen3)*
     &           (3*CW*dup2*USfC(Sfe4,1,tM1,Gen3) + 
     &             2*dup3*USfC(Sfe4,2,tM1,Gen3))))/
     &      (CW2*MW2*SB2*SW2)

	ENDLOOP(Sfe3)
	ENDLOOP(Sfe4)
	ENDLOOP(Gen3)

#ifdef DETAILED_DEBUG
	DHIGGS "dse_HHHH_mfv =", dse_HHHH_mfv ENDL
#endif

	if( mssmpart .le. 1 ) return

	LOOP(Gen3, gM,3,1)

        dse_HHHH_mfv = dse_HHHH_mfv - 
     &    3/(4.D0*Pi)*(Alfa1L*CA2*Mf2(bX,Gen3)*
     &        (B1(k2,Mf2(bX,Gen3),Mf2(bX,Gen3)) + 
     &          k2*DB1(k2,Mf2(bX,Gen3),Mf2(bX,Gen3)) + 
     &          2*DB0(k2,Mf2(bX,Gen3),Mf2(bX,Gen3))*Mf2(bX,Gen3)))/
     &      (CB2*MW2*SW2)

	ENDLOOP(Gen3)

	LOOP(Gen3, gM,3,1)
	LOOP(Sfe4, 1,2,1)
	LOOP(Sfe3, 1,2,1)

	dup5 = CA*KfC(Gen3,Gen3,bX) - MUE*SA*Mf(bX,Gen3)

	dup6 = CA*Kf(Gen3,Gen3,bX) - MUEC*SA*Mf(bX,Gen3)

	dup7 = CAB*CB*MW*MZ*SW2 - 3*CA*CW*Mf2(bX,Gen3)

	dup8 = -(CAB*CB*MW*MZ*(3 - 2*SW2)) + 6*CA*CW*Mf2(bX,Gen3)

        dse_HHHH_mfv = dse_HHHH_mfv + 
     &    1/(48.D0*Pi)*(Alfa1L*
     &        DB0(k2,MSf2(Sfe3,bX,Gen3),MSf2(Sfe4,bX,Gen3))*
     &        (USf(Sfe4,1,bX,Gen3)*
     &           (dup8*USfC(Sfe3,1,bX,Gen3) + 
     &             3*CW*dup5*USfC(Sfe3,2,bX,Gen3)) + 
     &          USf(Sfe4,2,bX,Gen3)*
     &           (3*CW*dup6*USfC(Sfe3,1,bX,Gen3) - 
     &             2*dup7*USfC(Sfe3,2,bX,Gen3)))*
     &        (USf(Sfe3,1,bX,Gen3)*
     &           (dup8*USfC(Sfe4,1,bX,Gen3) + 
     &             3*CW*dup5*USfC(Sfe4,2,bX,Gen3)) + 
     &          USf(Sfe3,2,bX,Gen3)*
     &           (3*CW*dup6*USfC(Sfe4,1,bX,Gen3) - 
     &             2*dup7*USfC(Sfe4,2,bX,Gen3))))/(CB2*CW2*MW2*SW2)

	ENDLOOP(Sfe3)
	ENDLOOP(Sfe4)
	ENDLOOP(Gen3)

#ifdef DETAILED_DEBUG
	DHIGGS "dse_HHHH_mfv =", dse_HHHH_mfv ENDL
#endif

	if( mssmpart .le. 2 ) return

	LOOP(Gen3, 1,3,1)

        dse_HHHH_mfv = dse_HHHH_mfv + 
     &    1/(16.D0*Pi)*(Alfa1L*
     &        (CAB**2*CB2*MW2*MZ2*
     &           DB0(k2,MSf2(1,1,Gen3),MSf2(1,1,Gen3)) - 
     &          4*CA2*CW2*Mf2(2,Gen3)*
     &           (B1(k2,Mf2(2,Gen3),Mf2(2,Gen3)) + 
     &             k2*DB1(k2,Mf2(2,Gen3),Mf2(2,Gen3)) + 
     &             2*DB0(k2,Mf2(2,Gen3),Mf2(2,Gen3))*Mf2(2,Gen3))))
     &       /(CB2*CW2*MW2*SW2)

	ENDLOOP(Gen3)

	LOOP(Sfe4, 1,2,1)
	LOOP(Sfe3, 1,2,1)
	LOOP(Gen3, 1,3,1)

	dup9 = CA*KfC(Gen3,Gen3,2) - MUE*SA*Mf(2,Gen3)

	dup10 = CA*Kf(Gen3,Gen3,2) - MUEC*SA*Mf(2,Gen3)

	dup11 = -(CAB*CB*MW*MZ*SW2) + CA*CW*Mf2(2,Gen3)

	dup12 = -(CAB*CB*MW*MZ*(1 - 2*SW2)) + 2*CA*CW*Mf2(2,Gen3)

        dse_HHHH_mfv = dse_HHHH_mfv + 
     &    1/(16.D0*Pi)*(Alfa1L*
     &        DB0(k2,MSf2(Sfe3,2,Gen3),MSf2(Sfe4,2,Gen3))*
     &        (USf(Sfe4,2,2,Gen3)*
     &           (CW*dup10*USfC(Sfe3,1,2,Gen3) + 
     &             2*dup11*USfC(Sfe3,2,2,Gen3)) + 
     &          USf(Sfe4,1,2,Gen3)*
     &           (dup12*USfC(Sfe3,1,2,Gen3) + 
     &             CW*dup9*USfC(Sfe3,2,2,Gen3)))*
     &        (USf(Sfe3,2,2,Gen3)*
     &           (CW*dup10*USfC(Sfe4,1,2,Gen3) + 
     &             2*dup11*USfC(Sfe4,2,2,Gen3)) + 
     &          USf(Sfe3,1,2,Gen3)*
     &           (dup12*USfC(Sfe4,1,2,Gen3) + 
     &             CW*dup9*USfC(Sfe4,2,2,Gen3))))/(CB2*CW2*MW2*SW2)

	ENDLOOP(Gen3)
	ENDLOOP(Sfe3)
	ENDLOOP(Sfe4)

#ifdef DETAILED_DEBUG
	DHIGGS "dse_HHHH_mfv =", dse_HHHH_mfv ENDL
#endif

	if( mssmpart .le. 3 ) return

        dse_HHHH_mfv = dse_HHHH_mfv + 
     &    1/(32.D0*Pi)*(Alfa1L*
     &        (MW2*((C2A*CAB - 2*S2A*SAB)**2*DB0(k2,Mh02,Mh02) + 
     &             CAB**2*(C2B**2*DB0(k2,MA02,MA02) + 
     &                9*C2A**2*DB0(k2,MHH2,MHH2)) + 
     &             2*((2*CAB*S2A + C2A*SAB)**2*DB0(k2,Mh02,MHH2) + 
     &                (4*CBA2*CW2**2 + 
     &                   C2B*(C2B*CAB**2 - 4*CAB*CBA*CW2))*
     &                 DB0(k2,MHp2,MHp2))) + 
     &          4*(CAB*MW2*(CAB*S2B**2 - 2*CW2*S2B*SBA) - 
     &             CW2**2*(k2 + MHp2 - MW2)*SBA2)*DB0(k2,MHp2,MW2)+
     &            (C2B**2*CAB**2*MW2 + 
     &             2*CBA2*(7*MW2 - CW2*(k2 + MZ2)))*DB0(k2,MZ2,MZ2)
     &            + 2*((CAB**2*MW2*S2B**2 - CW2*(k2 + MA02)*SBA2)*
     &              DB0(k2,MA02,MZ2) + 
     &             (-(CBA2*CW2**2*(2*k2 - 12*MW2)) + 
     &                C2B**2*CAB**2*MW2)*DB0(k2,MW2,MW2) - 
     &             CW2*(SBA2*
     &                 (B0(k2,MA02,MZ2) - 
     &                   2*
     &                    (B1(k2,MA02,MZ2) + k2*DB1(k2,MA02,MZ2))+
     &                     CW2*
     &                    (2*B0(k2,MHp2,MW2) - 
     &                      4*
     &                       (B1(k2,MHp2,MW2) + 
     &                       k2*DB1(k2,MHp2,MW2)))) + 
     &                CBA2*
     &                 (B0(k2,MZ2,MZ2) + 
     &                   CW2*
     &                    (2*B0(k2,MW2,MW2) - 
     &                      4*(B1(k2,MW2,MW2) + k2*DB1(k2,MW2,MW2))
     &                      ) - 
     &                   2*(B1(k2,MZ2,MZ2) + k2*DB1(k2,MZ2,MZ2)))))
     &          ))/(CW2**2*SW2)

	LOOP(Cha4, 1,2,1)
	LOOP(Cha3, 1,2,1)

        dup13 = CA*UCha(Cha4,2)*VCha(Cha3,1) + 
     &    SA*UCha(Cha4,1)*VCha(Cha3,2)

        dup14 = CA*UCha(Cha3,2)*VCha(Cha4,1) + 
     &    SA*UCha(Cha3,1)*VCha(Cha4,2)

        dup15 = CA*UChaC(Cha4,2)*VChaC(Cha3,1) + 
     &    SA*UChaC(Cha4,1)*VChaC(Cha3,2)

        dup16 = CA*UChaC(Cha3,2)*VChaC(Cha4,1) + 
     &    SA*UChaC(Cha3,1)*VChaC(Cha4,2)

        dse_HHHH_mfv = dse_HHHH_mfv - 
     &    1/(4.D0*Pi)*(Alfa1L*
     &        ((dup13*dup15 + dup14*dup16)*
     &           (B1(k2,MCha2(Cha3),MCha2(Cha4)) + 
     &             k2*DB1(k2,MCha2(Cha3),MCha2(Cha4))) + 
     &          DB0(k2,MCha2(Cha3),MCha2(Cha4))*
     &           (dup16*(dup15*MCha(Cha3)*MCha(Cha4) + 
     &                dup14*MCha2(Cha3)) + 
     &             dup13*(dup14*MCha(Cha3)*MCha(Cha4) + 
     &                dup15*MCha2(Cha3)))))/SW2

	ENDLOOP(Cha3)
	ENDLOOP(Cha4)

	LOOP(Neu4, 1,4,1)
	LOOP(Neu3, 1,4,1)

	dup17 = SW*ZNeu(Neu3,1) - CW*ZNeu(Neu3,2)

	dup18 = CA*ZNeu(Neu3,3) - SA*ZNeu(Neu3,4)

	dup19 = SW*ZNeu(Neu4,1) - CW*ZNeu(Neu4,2)

	dup20 = CA*ZNeu(Neu4,3) - SA*ZNeu(Neu4,4)

	dup21 = SW*ZNeuC(Neu3,1) - CW*ZNeuC(Neu3,2)

	dup22 = CA*ZNeuC(Neu3,3) - SA*ZNeuC(Neu3,4)

	dup23 = SW*ZNeuC(Neu4,1) - CW*ZNeuC(Neu4,2)

	dup24 = CA*ZNeuC(Neu4,3) - SA*ZNeuC(Neu4,4)

        dse_HHHH_mfv = dse_HHHH_mfv - 
     &    1/(16.D0*Pi)*(Alfa1L*
     &        (2*(dup18*dup19 + dup17*dup20)*
     &           (dup22*dup23 + dup21*dup24)*
     &           (B1(k2,MNeu2(Neu3),MNeu2(Neu4)) + 
     &             k2*DB1(k2,MNeu2(Neu3),MNeu2(Neu4))) + 
     &          DB0(k2,MNeu2(Neu3),MNeu2(Neu4))*
     &           ((dup17**2*dup20**2 + 
     &                dup18*
     &                 (dup18*dup19**2 + 2*dup17*dup19*dup20) + 
     &                dup21**2*dup24**2 + 
     &                dup22*(dup22*dup23**2 + 2*dup21*dup23*dup24))
     &               *MNeu(Neu3)*MNeu(Neu4) + 
     &             2*(dup18*dup19 + dup17*dup20)*
     &              (dup22*dup23 + dup21*dup24)*MNeu2(Neu3))))/
     &      (CW2*SW2)

	ENDLOOP(Neu3)
	ENDLOOP(Neu4)

#ifdef DETAILED_DEBUG
	DHIGGS "dse_HHHH_mfv =", dse_HHHH_mfv ENDL
#endif

	end


************************************************************************


	ComplexType function dse_HHHH_nmfv(k2,bX)
	implicit none
	integer bX
	RealType k2

#include "FH.h"
#include "looptools.h"

	integer All3, All4, Cha3, Cha4, Gen3, Ind1, Ind2, Ind3, Ind4
	integer Neu3, Neu4

	ComplexType dup1, dup10, dup11, dup12, dup2, dup3, dup4, dup5
	ComplexType dup6, dup7, dup8, dup9, tmp1, tmp2, tmp3

	dse_HHHH_nmfv = 0

	LOOP(Gen3, gM,3,1)

        dse_HHHH_nmfv = dse_HHHH_nmfv - 
     &    3/(4.D0*Pi)*(Alfa1L*SA2*Mf2(tM1,Gen3)*
     &        (B1(k2,Mf2(tM1,Gen3),Mf2(tM1,Gen3)) + 
     &          k2*DB1(k2,Mf2(tM1,Gen3),Mf2(tM1,Gen3)) + 
     &          2*DB0(k2,Mf2(tM1,Gen3),Mf2(tM1,Gen3))*Mf2(tM1,Gen3)
     &          ))/(MW2*SB2*SW2)

	ENDLOOP(Gen3)

	LOOP(All4, 1,6,gM)
	LOOP(All3, 1,6,gM)

	tmp1 = DB0(k2,MASf2(All3,tM1),MASf2(All4,tM1))

	LOOP(Ind4, 1,3,1)
	LOOP(Ind3, 1,3,1)
	LOOP(Ind2, 1,3,1)
	LOOP(Ind1, 1,3,1)

        dse_HHHH_nmfv = dse_HHHH_nmfv + 
     &    1/(48.D0*Pi)*(Alfa1L*tmp1*
     &        ((Delta(Ind3,Ind4)*
     &              (CAB*MW*MZ*SB*(3 - 4*SW2) + 
     &                6*CW*SA*Mf2(tM1,Ind3))*UASf(All4,Ind4,tM1) + 
     &             3*CW*(SA*Kf(Ind3,Ind4,tM1) - 
     &                CA*MUEC*Delta(Ind3,Ind4)*Mf(tM1,Ind3))*
     &              UASf(All4,3 + Ind4,tM1))*UASfC(All3,Ind3,tM1)+
     &            (3*CW*(SA*KfC(Ind4,Ind3,tM1) - 
     &                CA*MUE*Delta(Ind3,Ind4)*Mf(tM1,Ind3))*
     &              UASf(All4,Ind4,tM1) + 
     &             2*Delta(Ind3,Ind4)*
     &              (2*CAB*MW*MZ*SB*SW2 + 3*CW*SA*Mf2(tM1,Ind3))*
     &              UASf(All4,3 + Ind4,tM1))*
     &           UASfC(All3,3 + Ind3,tM1))*
     &        ((Delta(Ind1,Ind2)*
     &              (CAB*MW*MZ*SB*(3 - 4*SW2) + 
     &                6*CW*SA*Mf2(tM1,Ind1))*UASf(All3,Ind2,tM1) + 
     &             3*CW*(SA*Kf(Ind1,Ind2,tM1) - 
     &                CA*MUEC*Delta(Ind1,Ind2)*Mf(tM1,Ind1))*
     &              UASf(All3,3 + Ind2,tM1))*UASfC(All4,Ind1,tM1)+
     &            (3*CW*(SA*KfC(Ind2,Ind1,tM1) - 
     &                CA*MUE*Delta(Ind1,Ind2)*Mf(tM1,Ind1))*
     &              UASf(All3,Ind2,tM1) + 
     &             2*Delta(Ind1,Ind2)*
     &              (2*CAB*MW*MZ*SB*SW2 + 3*CW*SA*Mf2(tM1,Ind1))*
     &              UASf(All3,3 + Ind2,tM1))*
     &           UASfC(All4,3 + Ind1,tM1)))/(CW2*MW2*SB2*SW2)

	ENDLOOP(Ind1)
	ENDLOOP(Ind2)
	ENDLOOP(Ind3)
	ENDLOOP(Ind4)

	ENDLOOP(All3)
	ENDLOOP(All4)

#ifdef DETAILED_DEBUG
	DHIGGS "dse_HHHH_nmfv =", dse_HHHH_nmfv ENDL
#endif

	if( mssmpart .le. 1 ) return

	LOOP(Gen3, gM,3,1)

        dse_HHHH_nmfv = dse_HHHH_nmfv - 
     &    3/(4.D0*Pi)*(Alfa1L*CA2*Mf2(bX,Gen3)*
     &        (B1(k2,Mf2(bX,Gen3),Mf2(bX,Gen3)) + 
     &          k2*DB1(k2,Mf2(bX,Gen3),Mf2(bX,Gen3)) + 
     &          2*DB0(k2,Mf2(bX,Gen3),Mf2(bX,Gen3))*Mf2(bX,Gen3)))/
     &      (CB2*MW2*SW2)

	ENDLOOP(Gen3)

	LOOP(All4, 1,6,gM)
	LOOP(All3, 1,6,gM)

	tmp2 = DB0(k2,MASf2(All3,bX),MASf2(All4,bX))

	LOOP(Ind4, 1,3,1)
	LOOP(Ind3, 1,3,1)
	LOOP(Ind2, 1,3,1)
	LOOP(Ind1, 1,3,1)

        dse_HHHH_nmfv = dse_HHHH_nmfv + 
     &    1/(48.D0*Pi)*(Alfa1L*tmp2*
     &        ((Delta(Ind3,Ind4)*
     &              (CAB*CB*MW*MZ*(3 - 2*SW2) - 
     &                6*CA*CW*Mf2(bX,Ind3))*UASf(All4,Ind4,bX) - 
     &             3*CW*(CA*Kf(Ind3,Ind4,bX) - 
     &                MUEC*SA*Delta(Ind3,Ind4)*Mf(bX,Ind3))*
     &              UASf(All4,3 + Ind4,bX))*UASfC(All3,Ind3,bX) - 
     &          (3*CW*(CA*KfC(Ind4,Ind3,bX) - 
     &                MUE*SA*Delta(Ind3,Ind4)*Mf(bX,Ind3))*
     &              UASf(All4,Ind4,bX) - 
     &             2*Delta(Ind3,Ind4)*
     &              (CAB*CB*MW*MZ*SW2 - 3*CA*CW*Mf2(bX,Ind3))*
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
	DHIGGS "dse_HHHH_nmfv =", dse_HHHH_nmfv ENDL
#endif

	if( mssmpart .le. 2 ) return

	LOOP(Gen3, 1,3,1)

        dse_HHHH_nmfv = dse_HHHH_nmfv + 
     &    1/(16.D0*Pi)*(Alfa1L*
     &        (CAB**2*CB2*MW2*MZ2*
     &           DB0(k2,MASf2(Gen3,1),MASf2(Gen3,1)) - 
     &          4*CA2*CW2*Mf2(2,Gen3)*
     &           (B1(k2,Mf2(2,Gen3),Mf2(2,Gen3)) + 
     &             k2*DB1(k2,Mf2(2,Gen3),Mf2(2,Gen3)) + 
     &             2*DB0(k2,Mf2(2,Gen3),Mf2(2,Gen3))*Mf2(2,Gen3))))
     &       /(CB2*CW2*MW2*SW2)

	ENDLOOP(Gen3)

	LOOP(All4, 1,6,1)
	LOOP(All3, 1,6,1)

	tmp3 = DB0(k2,MASf2(All3,2),MASf2(All4,2))

	LOOP(Ind4, 1,3,1)
	LOOP(Ind3, 1,3,1)
	LOOP(Ind2, 1,3,1)
	LOOP(Ind1, 1,3,1)

        dse_HHHH_nmfv = dse_HHHH_nmfv + 
     &    1/(16.D0*Pi)*(Alfa1L*tmp3*
     &        ((Delta(Ind3,Ind4)*
     &              (CAB*CB*MW*MZ*(1 - 2*SW2) - 
     &                2*CA*CW*Mf2(2,Ind3))*UASf(All4,Ind4,2) - 
     &             CW*(CA*Kf(Ind3,Ind4,2) - 
     &                MUEC*SA*Delta(Ind3,Ind4)*Mf(2,Ind3))*
     &              UASf(All4,3 + Ind4,2))*UASfC(All3,Ind3,2) - 
     &          (CW*(CA*KfC(Ind4,Ind3,2) - 
     &                MUE*SA*Delta(Ind3,Ind4)*Mf(2,Ind3))*
     &              UASf(All4,Ind4,2) - 
     &             2*Delta(Ind3,Ind4)*
     &              (CAB*CB*MW*MZ*SW2 - CA*CW*Mf2(2,Ind3))*
     &              UASf(All4,3 + Ind4,2))*UASfC(All3,3 + Ind3,2))*
     &        ((Delta(Ind1,Ind2)*
     &              (CAB*CB*MW*MZ*(1 - 2*SW2) - 
     &                2*CA*CW*Mf2(2,Ind1))*UASf(All3,Ind2,2) - 
     &             CW*(CA*Kf(Ind1,Ind2,2) - 
     &                MUEC*SA*Delta(Ind1,Ind2)*Mf(2,Ind1))*
     &              UASf(All3,3 + Ind2,2))*UASfC(All4,Ind1,2) - 
     &          (CW*(CA*KfC(Ind2,Ind1,2) - 
     &                MUE*SA*Delta(Ind1,Ind2)*Mf(2,Ind1))*
     &              UASf(All3,Ind2,2) - 
     &             2*Delta(Ind1,Ind2)*
     &              (CAB*CB*MW*MZ*SW2 - CA*CW*Mf2(2,Ind1))*
     &              UASf(All3,3 + Ind2,2))*UASfC(All4,3 + Ind1,2)))
     &       /(CB2*CW2*MW2*SW2)

	ENDLOOP(Ind1)
	ENDLOOP(Ind2)
	ENDLOOP(Ind3)
	ENDLOOP(Ind4)

	ENDLOOP(All3)
	ENDLOOP(All4)

#ifdef DETAILED_DEBUG
	DHIGGS "dse_HHHH_nmfv =", dse_HHHH_nmfv ENDL
#endif

	if( mssmpart .le. 3 ) return

        dse_HHHH_nmfv = dse_HHHH_nmfv + 
     &    1/(32.D0*Pi)*(Alfa1L*
     &        (MW2*((C2A*CAB - 2*S2A*SAB)**2*DB0(k2,Mh02,Mh02) + 
     &             CAB**2*(C2B**2*DB0(k2,MA02,MA02) + 
     &                9*C2A**2*DB0(k2,MHH2,MHH2)) + 
     &             2*((2*CAB*S2A + C2A*SAB)**2*DB0(k2,Mh02,MHH2) + 
     &                (4*CBA2*CW2**2 + 
     &                   C2B*(C2B*CAB**2 - 4*CAB*CBA*CW2))*
     &                 DB0(k2,MHp2,MHp2))) + 
     &          4*(CAB*MW2*(CAB*S2B**2 - 2*CW2*S2B*SBA) - 
     &             CW2**2*(k2 + MHp2 - MW2)*SBA2)*DB0(k2,MHp2,MW2)+
     &            (C2B**2*CAB**2*MW2 + 
     &             2*CBA2*(7*MW2 - CW2*(k2 + MZ2)))*DB0(k2,MZ2,MZ2)
     &            + 2*((CAB**2*MW2*S2B**2 - CW2*(k2 + MA02)*SBA2)*
     &              DB0(k2,MA02,MZ2) + 
     &             (-(CBA2*CW2**2*(2*k2 - 12*MW2)) + 
     &                C2B**2*CAB**2*MW2)*DB0(k2,MW2,MW2) - 
     &             CW2*(SBA2*
     &                 (B0(k2,MA02,MZ2) - 
     &                   2*
     &                    (B1(k2,MA02,MZ2) + k2*DB1(k2,MA02,MZ2))+
     &                     CW2*
     &                    (2*B0(k2,MHp2,MW2) - 
     &                      4*
     &                       (B1(k2,MHp2,MW2) + 
     &                       k2*DB1(k2,MHp2,MW2)))) + 
     &                CBA2*
     &                 (B0(k2,MZ2,MZ2) + 
     &                   CW2*
     &                    (2*B0(k2,MW2,MW2) - 
     &                      4*(B1(k2,MW2,MW2) + k2*DB1(k2,MW2,MW2))
     &                      ) - 
     &                   2*(B1(k2,MZ2,MZ2) + k2*DB1(k2,MZ2,MZ2)))))
     &          ))/(CW2**2*SW2)

	LOOP(Cha4, 1,2,1)
	LOOP(Cha3, 1,2,1)

        dup1 = CA*UCha(Cha4,2)*VCha(Cha3,1) + 
     &    SA*UCha(Cha4,1)*VCha(Cha3,2)

        dup2 = CA*UCha(Cha3,2)*VCha(Cha4,1) + 
     &    SA*UCha(Cha3,1)*VCha(Cha4,2)

        dup3 = CA*UChaC(Cha4,2)*VChaC(Cha3,1) + 
     &    SA*UChaC(Cha4,1)*VChaC(Cha3,2)

        dup4 = CA*UChaC(Cha3,2)*VChaC(Cha4,1) + 
     &    SA*UChaC(Cha3,1)*VChaC(Cha4,2)

        dse_HHHH_nmfv = dse_HHHH_nmfv - 
     &    1/(4.D0*Pi)*(Alfa1L*
     &        ((dup1*dup3 + dup2*dup4)*
     &           (B1(k2,MCha2(Cha3),MCha2(Cha4)) + 
     &             k2*DB1(k2,MCha2(Cha3),MCha2(Cha4))) + 
     &          DB0(k2,MCha2(Cha3),MCha2(Cha4))*
     &           (dup4*(dup3*MCha(Cha3)*MCha(Cha4) + 
     &                dup2*MCha2(Cha3)) + 
     &             dup1*(dup2*MCha(Cha3)*MCha(Cha4) + 
     &                dup3*MCha2(Cha3)))))/SW2

	ENDLOOP(Cha3)
	ENDLOOP(Cha4)

	LOOP(Neu4, 1,4,1)
	LOOP(Neu3, 1,4,1)

	dup5 = SW*ZNeu(Neu3,1) - CW*ZNeu(Neu3,2)

	dup6 = CA*ZNeu(Neu3,3) - SA*ZNeu(Neu3,4)

	dup7 = SW*ZNeu(Neu4,1) - CW*ZNeu(Neu4,2)

	dup8 = CA*ZNeu(Neu4,3) - SA*ZNeu(Neu4,4)

	dup9 = SW*ZNeuC(Neu3,1) - CW*ZNeuC(Neu3,2)

	dup10 = CA*ZNeuC(Neu3,3) - SA*ZNeuC(Neu3,4)

	dup11 = SW*ZNeuC(Neu4,1) - CW*ZNeuC(Neu4,2)

	dup12 = CA*ZNeuC(Neu4,3) - SA*ZNeuC(Neu4,4)

        dse_HHHH_nmfv = dse_HHHH_nmfv - 
     &    1/(16.D0*Pi)*(Alfa1L*
     &        (2*(dup6*dup7 + dup5*dup8)*
     &           (dup10*dup11 + dup12*dup9)*
     &           (B1(k2,MNeu2(Neu3),MNeu2(Neu4)) + 
     &             k2*DB1(k2,MNeu2(Neu3),MNeu2(Neu4))) + 
     &          DB0(k2,MNeu2(Neu3),MNeu2(Neu4))*
     &           ((dup5**2*dup8**2 + 
     &                dup6*(dup6*dup7**2 + 2*dup5*dup7*dup8) + 
     &                dup12**2*dup9**2 + 
     &                dup10*(dup10*dup11**2 + 2*dup11*dup12*dup9))*
     &              MNeu(Neu3)*MNeu(Neu4) + 
     &             2*(dup6*dup7 + dup5*dup8)*
     &              (dup10*dup11 + dup12*dup9)*MNeu2(Neu3))))/
     &      (CW2*SW2)

	ENDLOOP(Neu3)
	ENDLOOP(Neu4)

#ifdef DETAILED_DEBUG
	DHIGGS "dse_HHHH_nmfv =", dse_HHHH_nmfv ENDL
#endif

	end

