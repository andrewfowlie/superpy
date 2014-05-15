* se_ZZ.F
* this file is part of FeynHiggs
* generated 05 Nov 2013 19:22

#include "externals.h"
#include "types.h"
#include "debug.h"


	ComplexType function se_ZZ(k2, bX)
	implicit none
	RealType k2
	integer bX


#include "FH.h"

	ComplexType se_ZZ_mfv, se_ZZ_nmfv
	external se_ZZ_mfv, se_ZZ_nmfv

	if( ibits(fv, 3, 2) .eq. 0 ) then
	  se_ZZ = se_ZZ_mfv(k2, bX)
	else
	  se_ZZ = se_ZZ_nmfv(k2, bX)
	endif

#ifdef DEBUG
	DHIGGS 'k2 =', k2 ENDL
	DHIGGS 'se_ZZ =', se_ZZ ENDL
#endif
	end


************************************************************************


	ComplexType function se_ZZ_mfv(k2,bX)
	implicit none
	integer bX
	RealType k2

#include "FH.h"
#include "looptools.h"

	integer Cha3, Cha4, Gen3, Neu3, Neu4, Sfe3, Sfe4

	ComplexType dup1, dup2, dup3, dup4, dup5, dup6, tmp1, tmp2

	se_ZZ_mfv = 0

	LOOP(Gen3, gM,3,1)

        se_ZZ_mfv = se_ZZ_mfv + 
     &    1/(24.D0*Pi)*(Alfa1L*
     &        ((9 - (24 - 32*SW2)*SW2)*
     &           (-A0(Mf2(tM1,Gen3)) + 
     &             2*B00(k2,Mf2(tM1,Gen3),Mf2(tM1,Gen3)) - 
     &             k2*B1(k2,Mf2(tM1,Gen3),Mf2(tM1,Gen3))) - 
     &          9*B0(k2,Mf2(tM1,Gen3),Mf2(tM1,Gen3))*Mf2(tM1,Gen3))
     &        )/(CW2*SW2)

	ENDLOOP(Gen3)

	LOOP(Gen3, gM,3,1)
	LOOP(Sfe3, 1,2,1)

        se_ZZ_mfv = se_ZZ_mfv + 
     &    1/(24.D0*Pi)*(Alfa1L*A0(MSf2(Sfe3,tM1,Gen3))*
     &        ((3 - 4*SW2)**2*USf2(Sfe3,1,tM1,Gen3) + 
     &          16*SW2**2*USf2(Sfe3,2,tM1,Gen3)))/(CW2*SW2)

	ENDLOOP(Sfe3)
	ENDLOOP(Gen3)

	LOOP(Gen3, gM,3,1)
	LOOP(Sfe4, 1,2,1)
	LOOP(Sfe3, 1,2,1)

        se_ZZ_mfv = se_ZZ_mfv - 
     &    1/(12.D0*Pi)*(Alfa1L*
     &        B00(k2,MSf2(Sfe3,tM1,Gen3),MSf2(Sfe4,tM1,Gen3))*
     &        (3*USf(Sfe4,1,tM1,Gen3)*USfC(Sfe3,1,tM1,Gen3) - 
     &          4*SW2*(USf(Sfe4,1,tM1,Gen3)*
     &              USfC(Sfe3,1,tM1,Gen3) + 
     &             USf(Sfe4,2,tM1,Gen3)*USfC(Sfe3,2,tM1,Gen3)))*
     &        (3*USf(Sfe3,1,tM1,Gen3)*USfC(Sfe4,1,tM1,Gen3) - 
     &          4*SW2*(USf(Sfe3,1,tM1,Gen3)*
     &              USfC(Sfe4,1,tM1,Gen3) + 
     &             USf(Sfe3,2,tM1,Gen3)*USfC(Sfe4,2,tM1,Gen3))))/
     &      (CW2*SW2)

	ENDLOOP(Sfe3)
	ENDLOOP(Sfe4)
	ENDLOOP(Gen3)

#ifdef DETAILED_DEBUG
	DHIGGS "se_ZZ_mfv =", se_ZZ_mfv ENDL
#endif

	if( mssmpart .le. 1 ) return

	LOOP(Gen3, gM,3,1)

        se_ZZ_mfv = se_ZZ_mfv + 
     &    1/(24.D0*Pi)*(Alfa1L*
     &        ((9 - (12 - 8*SW2)*SW2)*
     &           (-A0(Mf2(bX,Gen3)) + 
     &             2*B00(k2,Mf2(bX,Gen3),Mf2(bX,Gen3)) - 
     &             k2*B1(k2,Mf2(bX,Gen3),Mf2(bX,Gen3))) - 
     &          9*B0(k2,Mf2(bX,Gen3),Mf2(bX,Gen3))*Mf2(bX,Gen3)))/
     &      (CW2*SW2)

	ENDLOOP(Gen3)

	LOOP(Gen3, gM,3,1)
	LOOP(Sfe3, 1,2,1)

        se_ZZ_mfv = se_ZZ_mfv + 
     &    1/(24.D0*Pi)*(Alfa1L*A0(MSf2(Sfe3,bX,Gen3))*
     &        ((3 - 2*SW2)**2*USf2(Sfe3,1,bX,Gen3) + 
     &          4*SW2**2*USf2(Sfe3,2,bX,Gen3)))/(CW2*SW2)

	ENDLOOP(Sfe3)
	ENDLOOP(Gen3)

	LOOP(Gen3, gM,3,1)
	LOOP(Sfe4, 1,2,1)
	LOOP(Sfe3, 1,2,1)

        se_ZZ_mfv = se_ZZ_mfv - 
     &    1/(12.D0*Pi)*(Alfa1L*
     &        B00(k2,MSf2(Sfe3,bX,Gen3),MSf2(Sfe4,bX,Gen3))*
     &        (3*USf(Sfe4,1,bX,Gen3)*USfC(Sfe3,1,bX,Gen3) - 
     &          2*SW2*(USf(Sfe4,1,bX,Gen3)*USfC(Sfe3,1,bX,Gen3) + 
     &             USf(Sfe4,2,bX,Gen3)*USfC(Sfe3,2,bX,Gen3)))*
     &        (3*USf(Sfe3,1,bX,Gen3)*USfC(Sfe4,1,bX,Gen3) - 
     &          2*SW2*(USf(Sfe3,1,bX,Gen3)*USfC(Sfe4,1,bX,Gen3) + 
     &             USf(Sfe3,2,bX,Gen3)*USfC(Sfe4,2,bX,Gen3))))/
     &      (CW2*SW2)

	ENDLOOP(Sfe3)
	ENDLOOP(Sfe4)
	ENDLOOP(Gen3)

#ifdef DETAILED_DEBUG
	DHIGGS "se_ZZ_mfv =", se_ZZ_mfv ENDL
#endif

	if( mssmpart .le. 2 ) return

        se_ZZ_mfv = se_ZZ_mfv + 
     &    3/(8.D0*Pi)*(Alfa1L*(2*B00(k2,0.D0,0.D0) - k2*B1(k2,0.D0,0.D0)))/
     &      (CW2*SW2)

	LOOP(Gen3, 1,3,1)

        se_ZZ_mfv = se_ZZ_mfv + 
     &    1/(8.D0*Pi)*(Alfa1L*
     &        (A0(MSf2(1,1,Gen3)) - 
     &          2*B00(k2,MSf2(1,1,Gen3),MSf2(1,1,Gen3)) - 
     &          (1 - (4 - 8*SW2)*SW2)*
     &           (A0(Mf2(2,Gen3)) - 
     &             2*B00(k2,Mf2(2,Gen3),Mf2(2,Gen3)) + 
     &             k2*B1(k2,Mf2(2,Gen3),Mf2(2,Gen3))) - 
     &          B0(k2,Mf2(2,Gen3),Mf2(2,Gen3))*Mf2(2,Gen3)))/
     &      (CW2*SW2)

	ENDLOOP(Gen3)

	LOOP(Sfe3, 1,2,1)
	LOOP(Gen3, 1,3,1)

        se_ZZ_mfv = se_ZZ_mfv + 
     &    1/(8.D0*Pi)*(Alfa1L*A0(MSf2(Sfe3,2,Gen3))*
     &        ((1 - 2*SW2)**2*USf2(Sfe3,1,2,Gen3) + 
     &          4*SW2**2*USf2(Sfe3,2,2,Gen3)))/(CW2*SW2)

	ENDLOOP(Gen3)
	ENDLOOP(Sfe3)

	LOOP(Sfe4, 1,2,1)
	LOOP(Sfe3, 1,2,1)
	LOOP(Gen3, 1,3,1)

        se_ZZ_mfv = se_ZZ_mfv - 
     &    1/(4.D0*Pi)*(Alfa1L*
     &        B00(k2,MSf2(Sfe3,2,Gen3),MSf2(Sfe4,2,Gen3))*
     &        (USf(Sfe4,1,2,Gen3)*USfC(Sfe3,1,2,Gen3) - 
     &          2*SW2*(USf(Sfe4,1,2,Gen3)*USfC(Sfe3,1,2,Gen3) + 
     &             USf(Sfe4,2,2,Gen3)*USfC(Sfe3,2,2,Gen3)))*
     &        (USf(Sfe3,1,2,Gen3)*USfC(Sfe4,1,2,Gen3) - 
     &          2*SW2*(USf(Sfe3,1,2,Gen3)*USfC(Sfe4,1,2,Gen3) + 
     &             USf(Sfe3,2,2,Gen3)*USfC(Sfe4,2,2,Gen3))))/
     &      (CW2*SW2)

	ENDLOOP(Gen3)
	ENDLOOP(Sfe3)
	ENDLOOP(Sfe4)

#ifdef DETAILED_DEBUG
	DHIGGS "se_ZZ_mfv =", se_ZZ_mfv ENDL
#endif

	if( mssmpart .le. 3 ) return

        se_ZZ_mfv = se_ZZ_mfv + 
     &    1/(16.D0*Pi)*(Alfa1L*
     &        (4*MW2*(SBA2*B0(k2,Mh02,MZ2) + 
     &             CBA2*B0(k2,MHH2,MZ2)) + 
     &          CW2*(A0(MA02) + A0(Mh02) + A0(MHH2) + A0(MZ2) - 
     &             4*(SBA2*
     &                 (B00(k2,Mh02,MZ2) + B00(k2,MHH2,MA02)) + 
     &                CBA2*(B00(k2,Mh02,MA02) + B00(k2,MHH2,MZ2)))+
     &               2*((CW2*(9*CW2 - 2*SW2) + SW2**2)*A0(MW2) + 
     &                (CW2 - SW2)**2*
     &                 (A0(MHp2) - 2*B00(k2,MHp2,MHp2))) - 
     &             CW2*(((CW2**2*(20*k2 + 8*MW2) - 8*MW2*SW2**2)*
     &                    B0(k2,MW2,MW2) + 
     &                   4*(CW2*(9*CW2 - 2*SW2) + SW2**2)*
     &                    B00(k2,MW2,MW2))/CW2 + 
     &                8*CW2*k2*B1(k2,MW2,MW2)))))/(CW2**2*SW2)

	LOOP(Cha4, 1,2,1)

	tmp1 = A0(MCha2(Cha4))

	LOOP(Cha3, 1,2,1)

        dup1 = SW2*Delta(Cha3,Cha4) - 
     &    1/2.D0*(UCha(Cha4,2)*UChaC(Cha3,2)) - 
     &    UCha(Cha4,1)*UChaC(Cha3,1)

        dup2 = SW2*Delta(Cha3,Cha4) - 
     &    1/2.D0*(UCha(Cha3,2)*UChaC(Cha4,2)) - 
     &    UCha(Cha3,1)*UChaC(Cha4,1)

        dup3 = SW2*Delta(Cha3,Cha4) - 
     &    1/2.D0*(VCha(Cha4,2)*VChaC(Cha3,2)) - 
     &    VCha(Cha4,1)*VChaC(Cha3,1)

        dup4 = SW2*Delta(Cha3,Cha4) - 
     &    1/2.D0*(VCha(Cha3,2)*VChaC(Cha4,2)) - 
     &    VCha(Cha3,1)*VChaC(Cha4,1)

        se_ZZ_mfv = se_ZZ_mfv + 
     &    1/(8.D0*Pi)*(Alfa1L*
     &        (4*((dup1*dup2 + dup3*dup4)*
     &              (2*B00(k2,MCha2(Cha3),MCha2(Cha4)) - 
     &                k2*B1(k2,MCha2(Cha3),MCha2(Cha4))) + 
     &             B0(k2,MCha2(Cha3),MCha2(Cha4))*
     &              (dup1*(dup3*MCha(Cha3)*MCha(Cha4) - 
     &                   dup2*MCha2(Cha3)) + 
     &                dup4*
     &                 (dup2*MCha(Cha3)*MCha(Cha4) - 
     &                   dup3*MCha2(Cha3)))) - 
     &          tmp1*((2*SW2*Delta(Cha3,Cha4) - 
     &                2*UCha(Cha4,1)*UChaC(Cha3,1) - 
     &                UCha(Cha4,2)*UChaC(Cha3,2))*
     &              (2*SW2*Delta(Cha3,Cha4) - 
     &                2*UCha(Cha3,1)*UChaC(Cha4,1) - 
     &                UCha(Cha3,2)*UChaC(Cha4,2)) + 
     &             (2*SW2*Delta(Cha3,Cha4) - 
     &                2*VCha(Cha4,1)*VChaC(Cha3,1) - 
     &                VCha(Cha4,2)*VChaC(Cha3,2))*
     &              (2*SW2*Delta(Cha3,Cha4) - 
     &                2*VCha(Cha3,1)*VChaC(Cha4,1) - 
     &                VCha(Cha3,2)*VChaC(Cha4,2)))))/(CW2*SW2)

	ENDLOOP(Cha3)

	ENDLOOP(Cha4)

	LOOP(Neu4, 1,4,1)

	tmp2 = A0(MNeu2(Neu4))

	LOOP(Neu3, 1,4,1)

	dup5 = ZNeu(Neu4,3)*ZNeuC(Neu3,3) - ZNeu(Neu4,4)*ZNeuC(Neu3,4)

	dup6 = ZNeu(Neu3,3)*ZNeuC(Neu4,3) - ZNeu(Neu3,4)*ZNeuC(Neu4,4)

        se_ZZ_mfv = se_ZZ_mfv - 
     &    1/(16.D0*Pi)*(Alfa1L*
     &        (B0(k2,MNeu2(Neu3),MNeu2(Neu4))*
     &           (dup6**2*MNeu(Neu3)*MNeu(Neu4) + 
     &             dup5*(dup5*MNeu(Neu3)*MNeu(Neu4) + 
     &                2*dup6*MNeu2(Neu3))) - 
     &          dup6*(4*dup5*B00(k2,MNeu2(Neu3),MNeu2(Neu4)) - 
     &             2*tmp2*(ZNeu(Neu4,3)*ZNeuC(Neu3,3) - 
     &                ZNeu(Neu4,4)*ZNeuC(Neu3,4))) + 
     &          2*dup5*k2*B1(k2,MNeu2(Neu3),MNeu2(Neu4))*
     &           (ZNeu(Neu3,3)*ZNeuC(Neu4,3) - 
     &             ZNeu(Neu3,4)*ZNeuC(Neu4,4))))/(CW2*SW2)

	ENDLOOP(Neu3)

	ENDLOOP(Neu4)

#ifdef DETAILED_DEBUG
	DHIGGS "se_ZZ_mfv =", se_ZZ_mfv ENDL
#endif

	end


************************************************************************


	ComplexType function se_ZZ_nmfv(k2,bX)
	implicit none
	integer bX
	RealType k2

#include "FH.h"
#include "looptools.h"

	integer All3, All4, Cha3, Cha4, Gen3, Ind1, Ind2, Neu3, Neu4

	ComplexType dup1, dup2, dup3, dup4, dup5, dup6, tmp1, tmp2
	ComplexType tmp3, tmp4, tmp5, tmp6, tmp7, tmp8

	se_ZZ_nmfv = 0

	LOOP(Gen3, gM,3,1)

        se_ZZ_nmfv = se_ZZ_nmfv + 
     &    1/(24.D0*Pi)*(Alfa1L*
     &        ((9 - (24 - 32*SW2)*SW2)*
     &           (-A0(Mf2(tM1,Gen3)) + 
     &             2*B00(k2,Mf2(tM1,Gen3),Mf2(tM1,Gen3)) - 
     &             k2*B1(k2,Mf2(tM1,Gen3),Mf2(tM1,Gen3))) - 
     &          9*B0(k2,Mf2(tM1,Gen3),Mf2(tM1,Gen3))*Mf2(tM1,Gen3))
     &        )/(CW2*SW2)

	ENDLOOP(Gen3)

	LOOP(All3, 1,6,gM)

	tmp1 = A0(MASf2(All3,tM1))

	LOOP(Ind1, 1,3,1)

        se_ZZ_nmfv = se_ZZ_nmfv + 
     &    1/(24.D0*Pi)*(Alfa1L*tmp1*
     &        ((3 - 4*SW2)**2*UASf(All3,Ind1,tM1)*
     &           UASfC(All3,Ind1,tM1) + 
     &          16*SW2**2*UASf(All3,3 + Ind1,tM1)*
     &           UASfC(All3,3 + Ind1,tM1)))/(CW2*SW2)

	ENDLOOP(Ind1)

	ENDLOOP(All3)

	LOOP(All4, 1,6,gM)
	LOOP(All3, 1,6,gM)

	tmp2 = B00(k2,MASf2(All3,tM1),MASf2(All4,tM1))

	LOOP(Ind2, 1,3,1)
	LOOP(Ind1, 1,3,1)

        se_ZZ_nmfv = se_ZZ_nmfv - 
     &    1/(12.D0*Pi)*(Alfa1L*tmp2*
     &        (3*UASf(All4,Ind2,tM1)*UASfC(All3,Ind2,tM1) - 
     &          4*SW2*(UASf(All4,Ind2,tM1)*UASfC(All3,Ind2,tM1) + 
     &             UASf(All4,3 + Ind2,tM1)*UASfC(All3,3 + Ind2,tM1)
     &             ))*(3*UASf(All3,Ind1,tM1)*
     &           UASfC(All4,Ind1,tM1) - 
     &          4*SW2*(UASf(All3,Ind1,tM1)*UASfC(All4,Ind1,tM1) + 
     &             UASf(All3,3 + Ind1,tM1)*UASfC(All4,3 + Ind1,tM1)
     &             )))/(CW2*SW2)

	ENDLOOP(Ind1)
	ENDLOOP(Ind2)

	ENDLOOP(All3)
	ENDLOOP(All4)

#ifdef DETAILED_DEBUG
	DHIGGS "se_ZZ_nmfv =", se_ZZ_nmfv ENDL
#endif

	if( mssmpart .le. 1 ) return

	LOOP(Gen3, gM,3,1)

        se_ZZ_nmfv = se_ZZ_nmfv + 
     &    1/(24.D0*Pi)*(Alfa1L*
     &        ((9 - (12 - 8*SW2)*SW2)*
     &           (-A0(Mf2(bX,Gen3)) + 
     &             2*B00(k2,Mf2(bX,Gen3),Mf2(bX,Gen3)) - 
     &             k2*B1(k2,Mf2(bX,Gen3),Mf2(bX,Gen3))) - 
     &          9*B0(k2,Mf2(bX,Gen3),Mf2(bX,Gen3))*Mf2(bX,Gen3)))/
     &      (CW2*SW2)

	ENDLOOP(Gen3)

	LOOP(All3, 1,6,gM)

	tmp3 = A0(MASf2(All3,bX))

	LOOP(Ind1, 1,3,1)

        se_ZZ_nmfv = se_ZZ_nmfv + 
     &    1/(24.D0*Pi)*(Alfa1L*tmp3*
     &        ((3 - 2*SW2)**2*UASf(All3,Ind1,bX)*
     &           UASfC(All3,Ind1,bX) + 
     &          4*SW2**2*UASf(All3,3 + Ind1,bX)*
     &           UASfC(All3,3 + Ind1,bX)))/(CW2*SW2)

	ENDLOOP(Ind1)

	ENDLOOP(All3)

	LOOP(All4, 1,6,gM)
	LOOP(All3, 1,6,gM)

	tmp4 = B00(k2,MASf2(All3,bX),MASf2(All4,bX))

	LOOP(Ind2, 1,3,1)
	LOOP(Ind1, 1,3,1)

        se_ZZ_nmfv = se_ZZ_nmfv - 
     &    1/(12.D0*Pi)*(Alfa1L*tmp4*
     &        (3*UASf(All4,Ind2,bX)*UASfC(All3,Ind2,bX) - 
     &          2*SW2*(UASf(All4,Ind2,bX)*UASfC(All3,Ind2,bX) + 
     &             UASf(All4,3 + Ind2,bX)*UASfC(All3,3 + Ind2,bX)))
     &         *(3*UASf(All3,Ind1,bX)*UASfC(All4,Ind1,bX) - 
     &          2*SW2*(UASf(All3,Ind1,bX)*UASfC(All4,Ind1,bX) + 
     &             UASf(All3,3 + Ind1,bX)*UASfC(All4,3 + Ind1,bX)))
     &        )/(CW2*SW2)

	ENDLOOP(Ind1)
	ENDLOOP(Ind2)

	ENDLOOP(All3)
	ENDLOOP(All4)

#ifdef DETAILED_DEBUG
	DHIGGS "se_ZZ_nmfv =", se_ZZ_nmfv ENDL
#endif

	if( mssmpart .le. 2 ) return

        se_ZZ_nmfv = se_ZZ_nmfv + 
     &    3/(8.D0*Pi)*(Alfa1L*(2*B00(k2,0.D0,0.D0) - k2*B1(k2,0.D0,0.D0)))/
     &      (CW2*SW2)

	LOOP(Gen3, 1,3,1)

        se_ZZ_nmfv = se_ZZ_nmfv + 
     &    1/(8.D0*Pi)*(Alfa1L*
     &        (A0(MASf2(Gen3,1)) - 
     &          2*B00(k2,MASf2(Gen3,1),MASf2(Gen3,1)) - 
     &          (1 - (4 - 8*SW2)*SW2)*
     &           (A0(Mf2(2,Gen3)) - 
     &             2*B00(k2,Mf2(2,Gen3),Mf2(2,Gen3)) + 
     &             k2*B1(k2,Mf2(2,Gen3),Mf2(2,Gen3))) - 
     &          B0(k2,Mf2(2,Gen3),Mf2(2,Gen3))*Mf2(2,Gen3)))/
     &      (CW2*SW2)

	ENDLOOP(Gen3)

	LOOP(All3, 1,6,1)

	tmp5 = A0(MASf2(All3,2))

	LOOP(Ind1, 1,3,1)

        se_ZZ_nmfv = se_ZZ_nmfv + 
     &    1/(8.D0*Pi)*(Alfa1L*tmp5*
     &        ((1 - 2*SW2)**2*UASf(All3,Ind1,2)*
     &           UASfC(All3,Ind1,2) + 
     &          4*SW2**2*UASf(All3,3 + Ind1,2)*
     &           UASfC(All3,3 + Ind1,2)))/(CW2*SW2)

	ENDLOOP(Ind1)

	ENDLOOP(All3)

	LOOP(All4, 1,6,1)
	LOOP(All3, 1,6,1)

	tmp6 = B00(k2,MASf2(All3,2),MASf2(All4,2))

	LOOP(Ind2, 1,3,1)
	LOOP(Ind1, 1,3,1)

        se_ZZ_nmfv = se_ZZ_nmfv - 
     &    1/(4.D0*Pi)*(Alfa1L*tmp6*
     &        (UASf(All4,Ind2,2)*UASfC(All3,Ind2,2) - 
     &          2*SW2*(UASf(All4,Ind2,2)*UASfC(All3,Ind2,2) + 
     &             UASf(All4,3 + Ind2,2)*UASfC(All3,3 + Ind2,2)))*
     &        (UASf(All3,Ind1,2)*UASfC(All4,Ind1,2) - 
     &          2*SW2*(UASf(All3,Ind1,2)*UASfC(All4,Ind1,2) + 
     &             UASf(All3,3 + Ind1,2)*UASfC(All4,3 + Ind1,2))))/
     &      (CW2*SW2)

	ENDLOOP(Ind1)
	ENDLOOP(Ind2)

	ENDLOOP(All3)
	ENDLOOP(All4)

#ifdef DETAILED_DEBUG
	DHIGGS "se_ZZ_nmfv =", se_ZZ_nmfv ENDL
#endif

	if( mssmpart .le. 3 ) return

        se_ZZ_nmfv = se_ZZ_nmfv + 
     &    1/(16.D0*Pi)*(Alfa1L*
     &        (4*MW2*(SBA2*B0(k2,Mh02,MZ2) + 
     &             CBA2*B0(k2,MHH2,MZ2)) + 
     &          CW2*(A0(MA02) + A0(Mh02) + A0(MHH2) + A0(MZ2) - 
     &             4*(SBA2*
     &                 (B00(k2,Mh02,MZ2) + B00(k2,MHH2,MA02)) + 
     &                CBA2*(B00(k2,Mh02,MA02) + B00(k2,MHH2,MZ2)))+
     &               2*((CW2*(9*CW2 - 2*SW2) + SW2**2)*A0(MW2) + 
     &                (CW2 - SW2)**2*
     &                 (A0(MHp2) - 2*B00(k2,MHp2,MHp2))) - 
     &             CW2*(((CW2**2*(20*k2 + 8*MW2) - 8*MW2*SW2**2)*
     &                    B0(k2,MW2,MW2) + 
     &                   4*(CW2*(9*CW2 - 2*SW2) + SW2**2)*
     &                    B00(k2,MW2,MW2))/CW2 + 
     &                8*CW2*k2*B1(k2,MW2,MW2)))))/(CW2**2*SW2)

	LOOP(Cha4, 1,2,1)

	tmp7 = A0(MCha2(Cha4))

	LOOP(Cha3, 1,2,1)

        dup1 = SW2*Delta(Cha3,Cha4) - 
     &    1/2.D0*(UCha(Cha4,2)*UChaC(Cha3,2)) - 
     &    UCha(Cha4,1)*UChaC(Cha3,1)

        dup2 = SW2*Delta(Cha3,Cha4) - 
     &    1/2.D0*(UCha(Cha3,2)*UChaC(Cha4,2)) - 
     &    UCha(Cha3,1)*UChaC(Cha4,1)

        dup3 = SW2*Delta(Cha3,Cha4) - 
     &    1/2.D0*(VCha(Cha4,2)*VChaC(Cha3,2)) - 
     &    VCha(Cha4,1)*VChaC(Cha3,1)

        dup4 = SW2*Delta(Cha3,Cha4) - 
     &    1/2.D0*(VCha(Cha3,2)*VChaC(Cha4,2)) - 
     &    VCha(Cha3,1)*VChaC(Cha4,1)

        se_ZZ_nmfv = se_ZZ_nmfv + 
     &    1/(8.D0*Pi)*(Alfa1L*
     &        (4*((dup1*dup2 + dup3*dup4)*
     &              (2*B00(k2,MCha2(Cha3),MCha2(Cha4)) - 
     &                k2*B1(k2,MCha2(Cha3),MCha2(Cha4))) + 
     &             B0(k2,MCha2(Cha3),MCha2(Cha4))*
     &              (dup1*(dup3*MCha(Cha3)*MCha(Cha4) - 
     &                   dup2*MCha2(Cha3)) + 
     &                dup4*
     &                 (dup2*MCha(Cha3)*MCha(Cha4) - 
     &                   dup3*MCha2(Cha3)))) - 
     &          tmp7*((2*SW2*Delta(Cha3,Cha4) - 
     &                2*UCha(Cha4,1)*UChaC(Cha3,1) - 
     &                UCha(Cha4,2)*UChaC(Cha3,2))*
     &              (2*SW2*Delta(Cha3,Cha4) - 
     &                2*UCha(Cha3,1)*UChaC(Cha4,1) - 
     &                UCha(Cha3,2)*UChaC(Cha4,2)) + 
     &             (2*SW2*Delta(Cha3,Cha4) - 
     &                2*VCha(Cha4,1)*VChaC(Cha3,1) - 
     &                VCha(Cha4,2)*VChaC(Cha3,2))*
     &              (2*SW2*Delta(Cha3,Cha4) - 
     &                2*VCha(Cha3,1)*VChaC(Cha4,1) - 
     &                VCha(Cha3,2)*VChaC(Cha4,2)))))/(CW2*SW2)

	ENDLOOP(Cha3)

	ENDLOOP(Cha4)

	LOOP(Neu4, 1,4,1)

	tmp8 = A0(MNeu2(Neu4))

	LOOP(Neu3, 1,4,1)

	dup5 = ZNeu(Neu4,3)*ZNeuC(Neu3,3) - ZNeu(Neu4,4)*ZNeuC(Neu3,4)

	dup6 = ZNeu(Neu3,3)*ZNeuC(Neu4,3) - ZNeu(Neu3,4)*ZNeuC(Neu4,4)

        se_ZZ_nmfv = se_ZZ_nmfv - 
     &    1/(16.D0*Pi)*(Alfa1L*
     &        (B0(k2,MNeu2(Neu3),MNeu2(Neu4))*
     &           (dup6**2*MNeu(Neu3)*MNeu(Neu4) + 
     &             dup5*(dup5*MNeu(Neu3)*MNeu(Neu4) + 
     &                2*dup6*MNeu2(Neu3))) - 
     &          dup6*(4*dup5*B00(k2,MNeu2(Neu3),MNeu2(Neu4)) - 
     &             2*tmp8*(ZNeu(Neu4,3)*ZNeuC(Neu3,3) - 
     &                ZNeu(Neu4,4)*ZNeuC(Neu3,4))) + 
     &          2*dup5*k2*B1(k2,MNeu2(Neu3),MNeu2(Neu4))*
     &           (ZNeu(Neu3,3)*ZNeuC(Neu4,3) - 
     &             ZNeu(Neu3,4)*ZNeuC(Neu4,4))))/(CW2*SW2)

	ENDLOOP(Neu3)

	ENDLOOP(Neu4)

#ifdef DETAILED_DEBUG
	DHIGGS "se_ZZ_nmfv =", se_ZZ_nmfv ENDL
#endif

	end

