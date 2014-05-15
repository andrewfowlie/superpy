* dse_WW.F
* this file is part of FeynHiggs
* generated 05 Nov 2013 19:22

#include "externals.h"
#include "types.h"
#include "debug.h"


	ComplexType function dse_WW(k2, bX)
	implicit none
	RealType k2
	integer bX


#include "FH.h"

	ComplexType dse_WW_mfv, dse_WW_nmfv
	external dse_WW_mfv, dse_WW_nmfv

	if( ibits(fv, 3, 2) .eq. 0 ) then
	  dse_WW = dse_WW_mfv(k2, bX)
	else
	  dse_WW = dse_WW_nmfv(k2, bX)
	endif

#ifdef DEBUG
	DHIGGS 'k2 =', k2 ENDL
	DHIGGS 'dse_WW =', dse_WW ENDL
#endif
	end


************************************************************************


	ComplexType function dse_WW_mfv(k2,bX)
	implicit none
	integer bX
	RealType k2

#include "FH.h"
#include "looptools.h"

	integer Cha3, Gen3, Gen4, Neu3, Sfe3, Sfe4

	ComplexType dup1, dup2, dup3, dup4, tmp1, tmp10, tmp2, tmp3
	ComplexType tmp4, tmp5, tmp6, tmp7, tmp8, tmp9

	dse_WW_mfv = 0

	LOOP(Gen3, gM,3,1)

	tmp1 = B0(k2,0.D0,Mf2(tM1,Gen3))

	tmp2 = B1(k2,0.D0,Mf2(tM1,Gen3))

	tmp3 = DB0(k2,0.D0,Mf2(tM1,Gen3))

	tmp4 = DB00(k2,0.D0,Mf2(tM1,Gen3))

	tmp5 = DB1(k2,0.D0,Mf2(tM1,Gen3))

	LOOP(Gen4, gM,3,1)

        dse_WW_mfv = dse_WW_mfv + 
     &    3/(4.D0*Pi)*(Alfa1L*CKMin(Gen3,Gen4)*CKMinC(Gen3,Gen4)*
     &        (tmp1 + tmp2 + 2*tmp4 + k2*tmp5 + 
     &          tmp3*(k2 - Mf2(tM1,Gen3))))/SW2

	ENDLOOP(Gen4)

	ENDLOOP(Gen3)

	LOOP(Gen4, gM,3,1)
	LOOP(Gen3, gM,3,1)
	LOOP(Sfe4, 1,2,1)
	LOOP(Sfe3, 1,2,1)

        dse_WW_mfv = dse_WW_mfv - 
     &    3/(2.D0*Pi)*(Alfa1L*CKMin(Gen3,Gen4)*CKMinC(Gen3,Gen4)*
     &        DB00(k2,MSf2(Sfe3,tM1,Gen3),MSf2(Sfe4,bX,Gen4))*
     &        USf(Sfe4,1,bX,Gen4)*USf2(Sfe3,1,tM1,Gen3)*
     &        USfC(Sfe4,1,bX,Gen4))/SW2

	ENDLOOP(Sfe3)
	ENDLOOP(Sfe4)
	ENDLOOP(Gen3)
	ENDLOOP(Gen4)

#ifdef DETAILED_DEBUG
	DHIGGS "dse_WW_mfv =", dse_WW_mfv ENDL
#endif

	if( mssmpart .le. 1 ) return

	LOOP(Gen3, gM,3,1)

	tmp6 = B0(k2,0.D0,Mf2(tM1,Gen3))

	tmp7 = B1(k2,0.D0,Mf2(tM1,Gen3))

	tmp8 = DB0(k2,0.D0,Mf2(tM1,Gen3))

	tmp9 = DB00(k2,0.D0,Mf2(tM1,Gen3))

	tmp10 = DB1(k2,0.D0,Mf2(tM1,Gen3))

	LOOP(Gen4, gM,3,1)

        dse_WW_mfv = dse_WW_mfv - 
     &    3/(4.D0*Pi)*(Alfa1L*CKMin(Gen3,Gen4)*CKMinC(Gen3,Gen4)*
     &        (k2*tmp10 + tmp6 + tmp7 + 2*tmp9 - 
     &          B0(k2,Mf2(bX,Gen4),Mf2(tM1,Gen3)) - 
     &          B1(k2,Mf2(bX,Gen4),Mf2(tM1,Gen3)) - 
     &          2*DB00(k2,Mf2(bX,Gen4),Mf2(tM1,Gen3)) - 
     &          k2*DB1(k2,Mf2(bX,Gen4),Mf2(tM1,Gen3)) + 
     &          tmp8*(k2 - Mf2(tM1,Gen3)) - 
     &          DB0(k2,Mf2(bX,Gen4),Mf2(tM1,Gen3))*
     &           (k2 - Mf2(tM1,Gen3))))/SW2

	ENDLOOP(Gen4)

	ENDLOOP(Gen3)

#ifdef DETAILED_DEBUG
	DHIGGS "dse_WW_mfv =", dse_WW_mfv ENDL
#endif

	if( mssmpart .le. 2 ) return

	LOOP(Gen3, 1,3,1)

        dse_WW_mfv = dse_WW_mfv - 
     &    1/(4.D0*Pi)*(Alfa1L*
     &        (B1(k2,0.D0,Mf2(2,Gen3)) - 2*DB00(k2,0.D0,Mf2(2,Gen3)) + 
     &          k2*DB1(k2,0.D0,Mf2(2,Gen3))))/SW2

	ENDLOOP(Gen3)

	LOOP(Sfe3, 1,2,1)
	LOOP(Gen3, 1,3,1)

        dse_WW_mfv = dse_WW_mfv - 
     &    1/(2.D0*Pi)*(Alfa1L*
     &        DB00(k2,MSf2(1,1,Gen3),MSf2(Sfe3,2,Gen3))*
     &        USf2(Sfe3,1,2,Gen3))/SW2

	ENDLOOP(Gen3)
	ENDLOOP(Sfe3)

#ifdef DETAILED_DEBUG
	DHIGGS "dse_WW_mfv =", dse_WW_mfv ENDL
#endif

	if( mssmpart .le. 3 ) return

        dse_WW_mfv = dse_WW_mfv - 
     &    1/(4.D0*Pi)*(Alfa1L*
     &        (DB00(k2,MA02,MHp2) - 
     &          SBA2*(MW2*DB0(k2,Mh02,MW2) - DB00(k2,Mh02,MW2) - 
     &             DB00(k2,MHH2,MHp2)) - 
     &          CBA2*(MW2*DB0(k2,MHH2,MW2) - DB00(k2,Mh02,MHp2) - 
     &             DB00(k2,MHH2,MW2)) + DB00(k2,MZ2,MW2) + 
     &          SW2*(3*B0(k2,0.D0,MW2) + 
     &             (3*k2 - MW2)*DB0(k2,0.D0,MW2) + 
     &             ((CW2**2*(3*k2 + 2*MZ2) - MW2*SW2**2)*
     &                DB0(k2,MW2,MZ2))/(CW2*SW2) + 
     &             8*DB00(k2,MW2,0.D0) - 
     &             2*(B1(k2,MW2,0.D0) + k2*DB1(k2,MW2,0.D0))) + 
     &          CW2*(3*B0(k2,MW2,MZ2) - 
     &             2*(B1(k2,MW2,MZ2) - 4*DB00(k2,MW2,MZ2) + 
     &                k2*DB1(k2,MW2,MZ2)))))/SW2

	LOOP(Neu3, 1,4,1)
	LOOP(Cha3, 1,2,1)

        dup1 = 1/sqrt2*
     &    (sqrt2*(UChaC(Cha3,1)*ZNeu(Neu3,2)) + 
     &      UChaC(Cha3,2)*ZNeu(Neu3,3))

        dup2 = 1/sqrt2*
     &    (sqrt2*(VChaC(Cha3,1)*ZNeu(Neu3,2)) - 
     &      VChaC(Cha3,2)*ZNeu(Neu3,4))

        dup3 = 1/sqrt2*
     &    (sqrt2*(UCha(Cha3,1)*ZNeuC(Neu3,2)) + 
     &      UCha(Cha3,2)*ZNeuC(Neu3,3))

        dup4 = 1/sqrt2*
     &    (sqrt2*(VCha(Cha3,1)*ZNeuC(Neu3,2)) - 
     &      VCha(Cha3,2)*ZNeuC(Neu3,4))

        dse_WW_mfv = dse_WW_mfv + 
     &    1/(2.D0*Pi)*(Alfa1L*
     &        ((dup1*dup3 + dup2*dup4)*
     &           (B0(k2,MCha2(Cha3),MNeu2(Neu3)) + 
     &             B1(k2,MCha2(Cha3),MNeu2(Neu3)) + 
     &             2*DB00(k2,MCha2(Cha3),MNeu2(Neu3)) + 
     &             k2*DB1(k2,MCha2(Cha3),MNeu2(Neu3))) + 
     &          DB0(k2,MCha2(Cha3),MNeu2(Neu3))*
     &           ((dup1*dup2 + dup3*dup4)*MCha(Cha3)*MNeu(Neu3) + 
     &             (dup1*dup3 + dup2*dup4)*(k2 - MNeu2(Neu3)))))/
     &      SW2

	ENDLOOP(Cha3)
	ENDLOOP(Neu3)

#ifdef DETAILED_DEBUG
	DHIGGS "dse_WW_mfv =", dse_WW_mfv ENDL
#endif

	end


************************************************************************


	ComplexType function dse_WW_nmfv(k2,bX)
	implicit none
	integer bX
	RealType k2

#include "FH.h"
#include "looptools.h"

	integer All3, All4, Cha3, Gen3, Gen4, Ind1, Ind2, Ind3, Ind4
	integer Neu3

	ComplexType dup1, dup2, dup3, dup4, tmp1, tmp10, tmp11, tmp12
	ComplexType tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9

	dse_WW_nmfv = 0

	LOOP(Gen3, gM,3,1)

	tmp1 = B0(k2,0.D0,Mf2(tM1,Gen3))

	tmp2 = B1(k2,0.D0,Mf2(tM1,Gen3))

	tmp3 = DB0(k2,0.D0,Mf2(tM1,Gen3))

	tmp4 = DB00(k2,0.D0,Mf2(tM1,Gen3))

	tmp5 = DB1(k2,0.D0,Mf2(tM1,Gen3))

	LOOP(Gen4, gM,3,1)

        dse_WW_nmfv = dse_WW_nmfv + 
     &    3/(4.D0*Pi)*(Alfa1L*CKMin(Gen3,Gen4)*CKMinC(Gen3,Gen4)*
     &        (tmp1 + tmp2 + 2*tmp4 + k2*tmp5 + 
     &          tmp3*(k2 - Mf2(tM1,Gen3))))/SW2

	ENDLOOP(Gen4)

	ENDLOOP(Gen3)

	LOOP(All4, 1,6,gM)
	LOOP(All3, 1,6,gM)

	tmp6 = DB00(k2,MASf2(All3,tM1),MASf2(All4,bX))

	LOOP(Ind4, 1,3,1)
	LOOP(Ind3, 1,3,1)
	LOOP(Ind2, 1,3,1)
	LOOP(Ind1, 1,3,1)

        dse_WW_nmfv = dse_WW_nmfv - 
     &    3/(2.D0*Pi)*(Alfa1L*tmp6*CKMin(Ind1,Ind2)*
     &        CKMinC(Ind3,Ind4)*UASf(All3,Ind1,tM1)*
     &        UASf(All4,Ind4,bX)*UASfC(All3,Ind3,tM1)*
     &        UASfC(All4,Ind2,bX))/SW2

	ENDLOOP(Ind1)
	ENDLOOP(Ind2)
	ENDLOOP(Ind3)
	ENDLOOP(Ind4)

	ENDLOOP(All3)
	ENDLOOP(All4)

#ifdef DETAILED_DEBUG
	DHIGGS "dse_WW_nmfv =", dse_WW_nmfv ENDL
#endif

	if( mssmpart .le. 1 ) return

	LOOP(Gen3, gM,3,1)

	tmp7 = B0(k2,0.D0,Mf2(tM1,Gen3))

	tmp8 = B1(k2,0.D0,Mf2(tM1,Gen3))

	tmp9 = DB0(k2,0.D0,Mf2(tM1,Gen3))

	tmp10 = DB00(k2,0.D0,Mf2(tM1,Gen3))

	tmp11 = DB1(k2,0.D0,Mf2(tM1,Gen3))

	LOOP(Gen4, gM,3,1)

        dse_WW_nmfv = dse_WW_nmfv - 
     &    3/(4.D0*Pi)*(Alfa1L*CKMin(Gen3,Gen4)*CKMinC(Gen3,Gen4)*
     &        (2*tmp10 + k2*tmp11 + tmp7 + tmp8 - 
     &          B0(k2,Mf2(bX,Gen4),Mf2(tM1,Gen3)) - 
     &          B1(k2,Mf2(bX,Gen4),Mf2(tM1,Gen3)) - 
     &          2*DB00(k2,Mf2(bX,Gen4),Mf2(tM1,Gen3)) - 
     &          k2*DB1(k2,Mf2(bX,Gen4),Mf2(tM1,Gen3)) + 
     &          tmp9*(k2 - Mf2(tM1,Gen3)) - 
     &          DB0(k2,Mf2(bX,Gen4),Mf2(tM1,Gen3))*
     &           (k2 - Mf2(tM1,Gen3))))/SW2

	ENDLOOP(Gen4)

	ENDLOOP(Gen3)

#ifdef DETAILED_DEBUG
	DHIGGS "dse_WW_nmfv =", dse_WW_nmfv ENDL
#endif

	if( mssmpart .le. 2 ) return

	LOOP(Gen3, 1,3,1)

        dse_WW_nmfv = dse_WW_nmfv - 
     &    1/(4.D0*Pi)*(Alfa1L*
     &        (B1(k2,0.D0,Mf2(2,Gen3)) - 2*DB00(k2,0.D0,Mf2(2,Gen3)) + 
     &          k2*DB1(k2,0.D0,Mf2(2,Gen3))))/SW2

	ENDLOOP(Gen3)

	LOOP(Gen3, 1,3,1)
	LOOP(All3, 1,6,1)

	tmp12 = DB00(k2,MASf2(All3,2),MASf2(Gen3,1))

	LOOP(Ind2, 1,3,1)
	LOOP(Ind1, 1,3,1)

        dse_WW_nmfv = dse_WW_nmfv - 
     &    1/(2.D0*Pi)*(Alfa1L*tmp12*UASf(All3,Ind1,2)*
     &        UASf(Gen3,Ind2,1)*UASfC(All3,Ind2,2)*
     &        UASfC(Gen3,Ind1,1))/SW2

	ENDLOOP(Ind1)
	ENDLOOP(Ind2)

	ENDLOOP(All3)
	ENDLOOP(Gen3)

#ifdef DETAILED_DEBUG
	DHIGGS "dse_WW_nmfv =", dse_WW_nmfv ENDL
#endif

	if( mssmpart .le. 3 ) return

        dse_WW_nmfv = dse_WW_nmfv - 
     &    1/(4.D0*Pi)*(Alfa1L*
     &        ((CW2**2*(3*k2 + 2*MZ2) - MW2*SW2**2)*
     &           DB0(k2,MW2,MZ2) + 
     &          CW2*(DB00(k2,MA02,MHp2) - 
     &             SBA2*(MW2*DB0(k2,Mh02,MW2) - 
     &                DB00(k2,Mh02,MW2) - DB00(k2,MHH2,MHp2)) - 
     &             CBA2*(MW2*DB0(k2,MHH2,MW2) - 
     &                DB00(k2,Mh02,MHp2) - DB00(k2,MHH2,MW2)) + 
     &             DB00(k2,MZ2,MW2) + 
     &             SW2*(3*B0(k2,0.D0,MW2) + 
     &                (3*k2 - MW2)*DB0(k2,0.D0,MW2) + 
     &                8*DB00(k2,MW2,0.D0) - 
     &                2*(B1(k2,MW2,0.D0) + k2*DB1(k2,MW2,0.D0))) + 
     &             CW2*(3*B0(k2,MW2,MZ2) - 
     &                2*(B1(k2,MW2,MZ2) - 4*DB00(k2,MW2,MZ2) + 
     &                   k2*DB1(k2,MW2,MZ2))))))/(CW2*SW2)

	LOOP(Neu3, 1,4,1)
	LOOP(Cha3, 1,2,1)

        dup1 = 1/sqrt2*
     &    (sqrt2*(UChaC(Cha3,1)*ZNeu(Neu3,2)) + 
     &      UChaC(Cha3,2)*ZNeu(Neu3,3))

        dup2 = 1/sqrt2*
     &    (sqrt2*(VChaC(Cha3,1)*ZNeu(Neu3,2)) - 
     &      VChaC(Cha3,2)*ZNeu(Neu3,4))

        dup3 = 1/sqrt2*
     &    (sqrt2*(UCha(Cha3,1)*ZNeuC(Neu3,2)) + 
     &      UCha(Cha3,2)*ZNeuC(Neu3,3))

        dup4 = 1/sqrt2*
     &    (sqrt2*(VCha(Cha3,1)*ZNeuC(Neu3,2)) - 
     &      VCha(Cha3,2)*ZNeuC(Neu3,4))

        dse_WW_nmfv = dse_WW_nmfv + 
     &    1/(2.D0*Pi)*(Alfa1L*
     &        ((dup1*dup3 + dup2*dup4)*
     &           (B0(k2,MCha2(Cha3),MNeu2(Neu3)) + 
     &             B1(k2,MCha2(Cha3),MNeu2(Neu3)) + 
     &             2*DB00(k2,MCha2(Cha3),MNeu2(Neu3)) + 
     &             k2*DB1(k2,MCha2(Cha3),MNeu2(Neu3))) + 
     &          DB0(k2,MCha2(Cha3),MNeu2(Neu3))*
     &           ((dup1*dup2 + dup3*dup4)*MCha(Cha3)*MNeu(Neu3) + 
     &             (dup1*dup3 + dup2*dup4)*(k2 - MNeu2(Neu3)))))/
     &      SW2

	ENDLOOP(Cha3)
	ENDLOOP(Neu3)

#ifdef DETAILED_DEBUG
	DHIGGS "dse_WW_nmfv =", dse_WW_nmfv ENDL
#endif

	end

