* se_A0Z.F
* this file is part of FeynHiggs
* generated 05 Nov 2013 19:22

#include "externals.h"
#include "types.h"
#include "debug.h"


	ComplexType function se_A0Z(k2, bX)
	implicit none
	RealType k2
	integer bX


#include "FH.h"

	ComplexType se_A0Z_mfv, se_A0Z_nmfv
	external se_A0Z_mfv, se_A0Z_nmfv

	if( ibits(fv, 3, 2) .eq. 0 ) then
	  se_A0Z = se_A0Z_mfv(k2, bX)
	else
	  se_A0Z = se_A0Z_nmfv(k2, bX)
	endif

#ifdef DEBUG
	DHIGGS 'k2 =', k2 ENDL
	DHIGGS 'se_A0Z =', se_A0Z ENDL
#endif
	end


************************************************************************


	ComplexType function se_A0Z_mfv(k2,bX)
	implicit none
	integer bX
	RealType k2

#include "FH.h"
#include "looptools.h"

	integer Cha3, Cha4, Gen3, Neu3, Neu4, Sfe3, Sfe4

	ComplexType dup1, dup10, dup11, dup12, dup13, dup14, dup15
	ComplexType dup2, dup3, dup4, dup5, dup6, dup7, dup8, dup9

	se_A0Z_mfv = 0

	LOOP(Gen3, gM,3,1)

        se_A0Z_mfv = se_A0Z_mfv + 
     &    (3*cI)/(8.D0*Pi)*(Alfa1L*
     &        B0(k2,Mf2(tM1,Gen3),Mf2(tM1,Gen3))*Mf2(tM1,Gen3))/
     &      (CW*MW*SW2*TB)

	ENDLOOP(Gen3)

	LOOP(Gen3, gM,3,1)
	LOOP(Sfe4, 1,2,1)
	LOOP(Sfe3, 1,2,1)

        se_A0Z_mfv = se_A0Z_mfv + 
     &    cI/(16.D0*Pi)*(Alfa1L*
     &        (B0(k2,MSf2(Sfe3,tM1,Gen3),MSf2(Sfe4,tM1,Gen3)) + 
     &          2*B1(k2,MSf2(Sfe3,tM1,Gen3),MSf2(Sfe4,tM1,Gen3)))*
     &        ((Kf(Gen3,Gen3,tM1) + MUEC*TB*Mf(tM1,Gen3))*
     &           USf(Sfe4,2,tM1,Gen3)*USfC(Sfe3,1,tM1,Gen3) - 
     &          (KfC(Gen3,Gen3,tM1) + MUE*TB*Mf(tM1,Gen3))*
     &           USf(Sfe4,1,tM1,Gen3)*USfC(Sfe3,2,tM1,Gen3))*
     &        (3*USf(Sfe3,1,tM1,Gen3)*USfC(Sfe4,1,tM1,Gen3) - 
     &          4*SW2*(USf(Sfe3,1,tM1,Gen3)*
     &              USfC(Sfe4,1,tM1,Gen3) + 
     &             USf(Sfe3,2,tM1,Gen3)*USfC(Sfe4,2,tM1,Gen3))))/
     &      (CW*MW*SW2*TB)

	ENDLOOP(Sfe3)
	ENDLOOP(Sfe4)
	ENDLOOP(Gen3)

#ifdef DETAILED_DEBUG
	DHIGGS "se_A0Z_mfv =", se_A0Z_mfv ENDL
#endif

	if( mssmpart .le. 1 ) return

	LOOP(Gen3, gM,3,1)

        se_A0Z_mfv = se_A0Z_mfv - 
     &    (3*cI)/(8.D0*Pi)*(Alfa1L*TB*
     &        B0(k2,Mf2(bX,Gen3),Mf2(bX,Gen3))*Mf2(bX,Gen3))/
     &      (CW*MW*SW2)

	ENDLOOP(Gen3)

	LOOP(Gen3, gM,3,1)
	LOOP(Sfe4, 1,2,1)
	LOOP(Sfe3, 1,2,1)

        se_A0Z_mfv = se_A0Z_mfv - 
     &    cI/(16.D0*Pi)*(Alfa1L*
     &        (B0(k2,MSf2(Sfe3,bX,Gen3),MSf2(Sfe4,bX,Gen3)) + 
     &          2*B1(k2,MSf2(Sfe3,bX,Gen3),MSf2(Sfe4,bX,Gen3)))*
     &        ((TB*Kf(Gen3,Gen3,bX) + MUEC*Mf(bX,Gen3))*
     &           USf(Sfe4,2,bX,Gen3)*USfC(Sfe3,1,bX,Gen3) - 
     &          (TB*KfC(Gen3,Gen3,bX) + MUE*Mf(bX,Gen3))*
     &           USf(Sfe4,1,bX,Gen3)*USfC(Sfe3,2,bX,Gen3))*
     &        (3*USf(Sfe3,1,bX,Gen3)*USfC(Sfe4,1,bX,Gen3) - 
     &          2*SW2*(USf(Sfe3,1,bX,Gen3)*USfC(Sfe4,1,bX,Gen3) + 
     &             USf(Sfe3,2,bX,Gen3)*USfC(Sfe4,2,bX,Gen3))))/
     &      (CW*MW*SW2)

	ENDLOOP(Sfe3)
	ENDLOOP(Sfe4)
	ENDLOOP(Gen3)

#ifdef DETAILED_DEBUG
	DHIGGS "se_A0Z_mfv =", se_A0Z_mfv ENDL
#endif

	if( mssmpart .le. 2 ) return

	LOOP(Gen3, 1,3,1)

        se_A0Z_mfv = se_A0Z_mfv - 
     &    cI/(8.D0*Pi)*(Alfa1L*TB*B0(k2,Mf2(2,Gen3),Mf2(2,Gen3))*
     &        Mf2(2,Gen3))/(CW*MW*SW2)

	ENDLOOP(Gen3)

	LOOP(Sfe4, 1,2,1)
	LOOP(Sfe3, 1,2,1)
	LOOP(Gen3, 1,3,1)

        se_A0Z_mfv = se_A0Z_mfv - 
     &    cI/(16.D0*Pi)*(Alfa1L*
     &        (B0(k2,MSf2(Sfe3,2,Gen3),MSf2(Sfe4,2,Gen3)) + 
     &          2*B1(k2,MSf2(Sfe3,2,Gen3),MSf2(Sfe4,2,Gen3)))*
     &        ((TB*Kf(Gen3,Gen3,2) + MUEC*Mf(2,Gen3))*
     &           USf(Sfe4,2,2,Gen3)*USfC(Sfe3,1,2,Gen3) - 
     &          (TB*KfC(Gen3,Gen3,2) + MUE*Mf(2,Gen3))*
     &           USf(Sfe4,1,2,Gen3)*USfC(Sfe3,2,2,Gen3))*
     &        (USf(Sfe3,1,2,Gen3)*USfC(Sfe4,1,2,Gen3) - 
     &          2*SW2*(USf(Sfe3,1,2,Gen3)*USfC(Sfe4,1,2,Gen3) + 
     &             USf(Sfe3,2,2,Gen3)*USfC(Sfe4,2,2,Gen3))))/
     &      (CW*MW*SW2)

	ENDLOOP(Gen3)
	ENDLOOP(Sfe3)
	ENDLOOP(Sfe4)

#ifdef DETAILED_DEBUG
	DHIGGS "se_A0Z_mfv =", se_A0Z_mfv ENDL
#endif

	if( mssmpart .le. 3 ) return

        se_A0Z_mfv = se_A0Z_mfv - 
     &    cI/(16.D0*Pi)*(Alfa1L*MW*
     &        (C2B*(CBA*SAB*
     &              (B0(k2,MA02,Mh02) + 2*B1(k2,Mh02,MA02)) + 
     &             CAB*SBA*(B0(k2,MA02,MHH2) + 2*B1(k2,MHH2,MA02)))
     &            + SBA*((4*CBA + S2B*SAB)*B0(k2,Mh02,MZ2) + 
     &             2*(S2B*SAB*B1(k2,Mh02,MZ2) + 
     &                CBA*B1(k2,MZ2,Mh02))) - 
     &          CBA*((CAB*S2B + 4*SBA)*B0(k2,MHH2,MZ2) + 
     &             2*(CAB*S2B*B1(k2,MHH2,MZ2) + 
     &                SBA*B1(k2,MZ2,MHH2)))))/(CW*CW2*SW2)

	LOOP(Cha4, 1,2,1)
	LOOP(Cha3, 1,2,1)

        dup1 = SW2*Delta(Cha3,Cha4) - 
     &    1/2.D0*(UCha(Cha3,2)*UChaC(Cha4,2)) - 
     &    UCha(Cha3,1)*UChaC(Cha4,1)

        dup2 = SB*UCha(Cha4,2)*VCha(Cha3,1) + 
     &    CB*UCha(Cha4,1)*VCha(Cha3,2)

        dup3 = SW2*Delta(Cha3,Cha4) - 
     &    1/2.D0*(VCha(Cha4,2)*VChaC(Cha3,2)) - 
     &    VCha(Cha4,1)*VChaC(Cha3,1)

        dup4 = SB*UChaC(Cha3,2)*VChaC(Cha4,1) + 
     &    CB*UChaC(Cha3,1)*VChaC(Cha4,2)

        se_A0Z_mfv = se_A0Z_mfv - 
     &    cI/(2.D0*Pi*sqrt2)*
     &     (Alfa1L*((dup1*dup2 - dup3*dup4)*
     &           B0(k2,MCha2(Cha3),MCha2(Cha4))*MCha(Cha3) - 
     &          B1(k2,MCha2(Cha3),MCha2(Cha4))*
     &           (dup3*(dup4*MCha(Cha3) - dup2*MCha(Cha4)) - 
     &             dup1*(dup2*MCha(Cha3) - dup4*MCha(Cha4)))))/
     &      (CW*SW2)

	ENDLOOP(Cha3)
	ENDLOOP(Cha4)

	LOOP(Neu4, 1,4,1)
	LOOP(Neu3, 1,4,1)

	dup5 = SW*ZNeu(Neu3,1) - CW*ZNeu(Neu3,2)

	dup6 = SB*ZNeu(Neu3,3) - CB*ZNeu(Neu3,4)

	dup7 = SW*ZNeu(Neu4,1) - CW*ZNeu(Neu4,2)

	dup8 = SB*ZNeu(Neu4,3) - CB*ZNeu(Neu4,4)

	dup9 = SW*ZNeuC(Neu3,1) - CW*ZNeuC(Neu3,2)

	dup10 = SB*ZNeuC(Neu3,3) - CB*ZNeuC(Neu3,4)

        dup11 = ZNeu(Neu4,3)*ZNeuC(Neu3,3) - 
     &    ZNeu(Neu4,4)*ZNeuC(Neu3,4)

	dup12 = SW*ZNeuC(Neu4,1) - CW*ZNeuC(Neu4,2)

	dup13 = SB*ZNeuC(Neu4,3) - CB*ZNeuC(Neu4,4)

        dup14 = ZNeu(Neu3,3)*ZNeuC(Neu4,3) - 
     &    ZNeu(Neu3,4)*ZNeuC(Neu4,4)

        dup15 = dup10*dup11*dup12 + 
     &    dup14*(dup6*dup7 + dup5*dup8) + dup11*dup13*dup9

        se_A0Z_mfv = se_A0Z_mfv + 
     &    cI/(16.D0*Pi)*(Alfa1L*
     &        (dup15*B0(k2,MNeu2(Neu3),MNeu2(Neu4))*MNeu(Neu3) + 
     &          B1(k2,MNeu2(Neu3),MNeu2(Neu4))*
     &           (dup15*MNeu(Neu3) - 
     &             (dup10*dup12*dup14 + 
     &                dup11*(dup6*dup7 + dup5*dup8) + 
     &                dup13*dup14*dup9)*MNeu(Neu4))))/(CW2*SW2)

	ENDLOOP(Neu3)
	ENDLOOP(Neu4)

#ifdef DETAILED_DEBUG
	DHIGGS "se_A0Z_mfv =", se_A0Z_mfv ENDL
#endif

	end


************************************************************************


	ComplexType function se_A0Z_nmfv(k2,bX)
	implicit none
	integer bX
	RealType k2

#include "FH.h"
#include "looptools.h"

	integer All3, All4, Cha3, Cha4, Gen3, Ind1, Ind2, Ind3, Neu3
	integer Neu4

	ComplexType dup1, dup10, dup11, dup12, dup13, dup14, dup15
	ComplexType dup2, dup3, dup4, dup5, dup6, dup7, dup8, dup9
	ComplexType tmp1, tmp2, tmp3, tmp4, tmp5, tmp6

	se_A0Z_nmfv = 0

	LOOP(Gen3, gM,3,1)

        se_A0Z_nmfv = se_A0Z_nmfv + 
     &    (3*cI)/(8.D0*Pi)*(Alfa1L*
     &        B0(k2,Mf2(tM1,Gen3),Mf2(tM1,Gen3))*Mf2(tM1,Gen3))/
     &      (CW*MW*SW2*TB)

	ENDLOOP(Gen3)

	LOOP(All4, 1,6,gM)
	LOOP(All3, 1,6,gM)

	tmp1 = B0(k2,MASf2(All3,tM1),MASf2(All4,tM1))

	tmp2 = B1(k2,MASf2(All3,tM1),MASf2(All4,tM1))

	LOOP(Ind3, 1,3,1)
	LOOP(Ind2, 1,3,1)
	LOOP(Ind1, 1,3,1)

        se_A0Z_nmfv = se_A0Z_nmfv + 
     &    cI/(16.D0*Pi)*(Alfa1L*(tmp1 + 2*tmp2)*
     &        ((Kf(Ind2,Ind3,tM1) + 
     &             MUEC*TB*Delta(Ind2,Ind3)*Mf(tM1,Ind2))*
     &           UASf(All4,3 + Ind3,tM1)*UASfC(All3,Ind2,tM1) - 
     &          (KfC(Ind3,Ind2,tM1) + 
     &             MUE*TB*Delta(Ind2,Ind3)*Mf(tM1,Ind2))*
     &           UASf(All4,Ind3,tM1)*UASfC(All3,3 + Ind2,tM1))*
     &        (3*UASf(All3,Ind1,tM1)*UASfC(All4,Ind1,tM1) - 
     &          4*SW2*(UASf(All3,Ind1,tM1)*UASfC(All4,Ind1,tM1) + 
     &             UASf(All3,3 + Ind1,tM1)*UASfC(All4,3 + Ind1,tM1)
     &             )))/(CW*MW*SW2*TB)

	ENDLOOP(Ind1)
	ENDLOOP(Ind2)
	ENDLOOP(Ind3)

	ENDLOOP(All3)
	ENDLOOP(All4)

#ifdef DETAILED_DEBUG
	DHIGGS "se_A0Z_nmfv =", se_A0Z_nmfv ENDL
#endif

	if( mssmpart .le. 1 ) return

	LOOP(Gen3, gM,3,1)

        se_A0Z_nmfv = se_A0Z_nmfv - 
     &    (3*cI)/(8.D0*Pi)*(Alfa1L*TB*
     &        B0(k2,Mf2(bX,Gen3),Mf2(bX,Gen3))*Mf2(bX,Gen3))/
     &      (CW*MW*SW2)

	ENDLOOP(Gen3)

	LOOP(All4, 1,6,gM)
	LOOP(All3, 1,6,gM)

	tmp3 = B0(k2,MASf2(All3,bX),MASf2(All4,bX))

	tmp4 = B1(k2,MASf2(All3,bX),MASf2(All4,bX))

	LOOP(Ind3, 1,3,1)
	LOOP(Ind2, 1,3,1)
	LOOP(Ind1, 1,3,1)

        se_A0Z_nmfv = se_A0Z_nmfv - 
     &    cI/(16.D0*Pi)*(Alfa1L*(tmp3 + 2*tmp4)*
     &        ((TB*Kf(Ind2,Ind3,bX) + 
     &             MUEC*Delta(Ind2,Ind3)*Mf(bX,Ind2))*
     &           UASf(All4,3 + Ind3,bX)*UASfC(All3,Ind2,bX) - 
     &          (TB*KfC(Ind3,Ind2,bX) + 
     &             MUE*Delta(Ind2,Ind3)*Mf(bX,Ind2))*
     &           UASf(All4,Ind3,bX)*UASfC(All3,3 + Ind2,bX))*
     &        (3*UASf(All3,Ind1,bX)*UASfC(All4,Ind1,bX) - 
     &          2*SW2*(UASf(All3,Ind1,bX)*UASfC(All4,Ind1,bX) + 
     &             UASf(All3,3 + Ind1,bX)*UASfC(All4,3 + Ind1,bX)))
     &        )/(CW*MW*SW2)

	ENDLOOP(Ind1)
	ENDLOOP(Ind2)
	ENDLOOP(Ind3)

	ENDLOOP(All3)
	ENDLOOP(All4)

#ifdef DETAILED_DEBUG
	DHIGGS "se_A0Z_nmfv =", se_A0Z_nmfv ENDL
#endif

	if( mssmpart .le. 2 ) return

	LOOP(Gen3, 1,3,1)

        se_A0Z_nmfv = se_A0Z_nmfv - 
     &    cI/(8.D0*Pi)*(Alfa1L*TB*B0(k2,Mf2(2,Gen3),Mf2(2,Gen3))*
     &        Mf2(2,Gen3))/(CW*MW*SW2)

	ENDLOOP(Gen3)

	LOOP(All4, 1,6,1)
	LOOP(All3, 1,6,1)

	tmp5 = B0(k2,MASf2(All3,2),MASf2(All4,2))

	tmp6 = B1(k2,MASf2(All3,2),MASf2(All4,2))

	LOOP(Ind3, 1,3,1)
	LOOP(Ind2, 1,3,1)
	LOOP(Ind1, 1,3,1)

        se_A0Z_nmfv = se_A0Z_nmfv - 
     &    cI/(16.D0*Pi)*(Alfa1L*(tmp5 + 2*tmp6)*
     &        ((TB*Kf(Ind2,Ind3,2) + 
     &             MUEC*Delta(Ind2,Ind3)*Mf(2,Ind2))*
     &           UASf(All4,3 + Ind3,2)*UASfC(All3,Ind2,2) - 
     &          (TB*KfC(Ind3,Ind2,2) + 
     &             MUE*Delta(Ind2,Ind3)*Mf(2,Ind2))*
     &           UASf(All4,Ind3,2)*UASfC(All3,3 + Ind2,2))*
     &        (UASf(All3,Ind1,2)*UASfC(All4,Ind1,2) - 
     &          2*SW2*(UASf(All3,Ind1,2)*UASfC(All4,Ind1,2) + 
     &             UASf(All3,3 + Ind1,2)*UASfC(All4,3 + Ind1,2))))/
     &      (CW*MW*SW2)

	ENDLOOP(Ind1)
	ENDLOOP(Ind2)
	ENDLOOP(Ind3)

	ENDLOOP(All3)
	ENDLOOP(All4)

#ifdef DETAILED_DEBUG
	DHIGGS "se_A0Z_nmfv =", se_A0Z_nmfv ENDL
#endif

	if( mssmpart .le. 3 ) return

        se_A0Z_nmfv = se_A0Z_nmfv - 
     &    cI/(16.D0*Pi)*(Alfa1L*MW*
     &        (C2B*(CBA*SAB*
     &              (B0(k2,MA02,Mh02) + 2*B1(k2,Mh02,MA02)) + 
     &             CAB*SBA*(B0(k2,MA02,MHH2) + 2*B1(k2,MHH2,MA02)))
     &            + SBA*((4*CBA + S2B*SAB)*B0(k2,Mh02,MZ2) + 
     &             2*(S2B*SAB*B1(k2,Mh02,MZ2) + 
     &                CBA*B1(k2,MZ2,Mh02))) - 
     &          CBA*((CAB*S2B + 4*SBA)*B0(k2,MHH2,MZ2) + 
     &             2*(CAB*S2B*B1(k2,MHH2,MZ2) + 
     &                SBA*B1(k2,MZ2,MHH2)))))/(CW*CW2*SW2)

	LOOP(Cha4, 1,2,1)
	LOOP(Cha3, 1,2,1)

        dup1 = SW2*Delta(Cha3,Cha4) - 
     &    1/2.D0*(UCha(Cha3,2)*UChaC(Cha4,2)) - 
     &    UCha(Cha3,1)*UChaC(Cha4,1)

        dup2 = SB*UCha(Cha4,2)*VCha(Cha3,1) + 
     &    CB*UCha(Cha4,1)*VCha(Cha3,2)

        dup3 = SW2*Delta(Cha3,Cha4) - 
     &    1/2.D0*(VCha(Cha4,2)*VChaC(Cha3,2)) - 
     &    VCha(Cha4,1)*VChaC(Cha3,1)

        dup4 = SB*UChaC(Cha3,2)*VChaC(Cha4,1) + 
     &    CB*UChaC(Cha3,1)*VChaC(Cha4,2)

        se_A0Z_nmfv = se_A0Z_nmfv - 
     &    cI/(2.D0*Pi*sqrt2)*
     &     (Alfa1L*((dup1*dup2 - dup3*dup4)*
     &           B0(k2,MCha2(Cha3),MCha2(Cha4))*MCha(Cha3) - 
     &          B1(k2,MCha2(Cha3),MCha2(Cha4))*
     &           (dup3*(dup4*MCha(Cha3) - dup2*MCha(Cha4)) - 
     &             dup1*(dup2*MCha(Cha3) - dup4*MCha(Cha4)))))/
     &      (CW*SW2)

	ENDLOOP(Cha3)
	ENDLOOP(Cha4)

	LOOP(Neu4, 1,4,1)
	LOOP(Neu3, 1,4,1)

	dup5 = SW*ZNeu(Neu3,1) - CW*ZNeu(Neu3,2)

	dup6 = SB*ZNeu(Neu3,3) - CB*ZNeu(Neu3,4)

	dup7 = SW*ZNeu(Neu4,1) - CW*ZNeu(Neu4,2)

	dup8 = SB*ZNeu(Neu4,3) - CB*ZNeu(Neu4,4)

	dup9 = SW*ZNeuC(Neu3,1) - CW*ZNeuC(Neu3,2)

	dup10 = SB*ZNeuC(Neu3,3) - CB*ZNeuC(Neu3,4)

        dup11 = ZNeu(Neu4,3)*ZNeuC(Neu3,3) - 
     &    ZNeu(Neu4,4)*ZNeuC(Neu3,4)

	dup12 = SW*ZNeuC(Neu4,1) - CW*ZNeuC(Neu4,2)

	dup13 = SB*ZNeuC(Neu4,3) - CB*ZNeuC(Neu4,4)

        dup14 = ZNeu(Neu3,3)*ZNeuC(Neu4,3) - 
     &    ZNeu(Neu3,4)*ZNeuC(Neu4,4)

        dup15 = dup10*dup11*dup12 + 
     &    dup14*(dup6*dup7 + dup5*dup8) + dup11*dup13*dup9

        se_A0Z_nmfv = se_A0Z_nmfv + 
     &    cI/(16.D0*Pi)*(Alfa1L*
     &        (dup15*B0(k2,MNeu2(Neu3),MNeu2(Neu4))*MNeu(Neu3) + 
     &          B1(k2,MNeu2(Neu3),MNeu2(Neu4))*
     &           (dup15*MNeu(Neu3) - 
     &             (dup10*dup12*dup14 + 
     &                dup11*(dup6*dup7 + dup5*dup8) + 
     &                dup13*dup14*dup9)*MNeu(Neu4))))/(CW2*SW2)

	ENDLOOP(Neu3)
	ENDLOOP(Neu4)

#ifdef DETAILED_DEBUG
	DHIGGS "se_A0Z_nmfv =", se_A0Z_nmfv ENDL
#endif

	end

