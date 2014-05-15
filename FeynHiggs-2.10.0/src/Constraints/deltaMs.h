
        cVLLSM = 1/Pi**2*
     &    (GF**2*CKM(3,3)**2*CKMC(3,2)**2*
     &      (-((D00z(MW2,MW2,Mf2(tT,3),Mf2(tT,3)) - 
     &             2*MW2*D0z(MW2,MW2,Mf2(tT,3),Mf2(tT,3)))*
     &           Mf2(tT,3)**2) + 
     &        MW2**2*(-2*(C0z(0.D0,0.D0,MW2) - 
     &              2*C0z(0.D0,MW2,Mf2(tT,3)) + 
     &              C0z(MW2,Mf2(tT,3),Mf2(tT,3)) + 
     &              12*D00z(0.D0,MW2,MW2,Mf2(tT,3)) - 
     &              6*(D00z(0.D0,0.D0,MW2,MW2) + 
     &                 D00z(MW2,MW2,Mf2(tT,3),Mf2(tT,3))) + 
     &              MW2*(D0z(0.D0,0.D0,MW2,MW2) - 
     &                 2*D0z(0.D0,MW2,MW2,Mf2(tT,3)) + 
     &                 D0z(MW2,MW2,Mf2(tT,3),Mf2(tT,3)))) + 
     &           2*(D0z(0.D0,MW2,MW2,Mf2(tT,3)) - 
     &              D0z(MW2,MW2,Mf2(tT,3),Mf2(tT,3)))*Mf2(tT,3))))

#ifdef DETAILED_DEBUG
	DCONST "cVLLSM =", cVLLSM ENDL
#endif


        cVLLHp = -(1/Pi**2*
     &      (GF**2*CKM(3,3)**2*CKMC(3,2)**2*
     &         (D00z(MHp2,MHp2,Mf2(tT,3),Mf2(tT,3)) + 
     &           2*TB2*(D00z(MHp2,MW2,Mf2(tT,3),Mf2(tT,3)) - 
     &              MW2*D0z(MHp2,MW2,Mf2(tT,3),Mf2(tT,3))))*
     &         Mf2(tT,3)**2)/TB2**2)

#ifdef DETAILED_DEBUG
	DCONST "cVLLHp =", cVLLHp ENDL
#endif


	cVLLCha = 0

	LOOP(Cha6, 1,2,1)
	LOOP(Cha5, 1,2,1)
	LOOP(All6, 1,6,1)
	LOOP(All5, 1,6,1)

        tmp1 = D00z(MASf2(All5,3),MASf2(All6,3),MCha2(Cha5),
     &    MCha2(Cha6))

	LOOP(Ind4, 1,3,1)
	LOOP(Ind3, 1,3,1)
	LOOP(Ind2, 1,3,1)
	LOOP(Ind1, 1,3,1)

        cVLLCha = cVLLCha - 
     &    1/(4.D0*Pi**2)*(GF**2*tmp1*CKM(Ind1,3)*CKM(Ind2,3)*
     &        CKMC(Ind3,2)*CKMC(Ind4,2)*
     &        (sqrt2*(Mf(3,Ind3)*UASfC(All5,3 + Ind3,3)*
     &             VCha(Cha5,2)) - 
     &          2*MW*SB*UASfC(All5,Ind3,3)*VCha(Cha5,1))*
     &        (sqrt2*(Mf(3,Ind4)*UASfC(All6,3 + Ind4,3)*
     &             VCha(Cha6,2)) - 
     &          2*MW*SB*UASfC(All6,Ind4,3)*VCha(Cha6,1))*
     &        (sqrt2*(Mf(3,Ind2)*UASf(All6,3 + Ind2,3)*
     &             VChaC(Cha5,2)) - 
     &          2*MW*SB*UASf(All6,Ind2,3)*VChaC(Cha5,1))*
     &        (sqrt2*(Mf(3,Ind1)*UASf(All5,3 + Ind1,3)*
     &             VChaC(Cha6,2)) - 
     &          2*MW*SB*UASf(All5,Ind1,3)*VChaC(Cha6,1)))/SB2**2

	ENDLOOP(Ind1)
	ENDLOOP(Ind2)
	ENDLOOP(Ind3)
	ENDLOOP(Ind4)

	ENDLOOP(All5)
	ENDLOOP(All6)
	ENDLOOP(Cha5)
	ENDLOOP(Cha6)

#ifdef DETAILED_DEBUG
	DCONST "cVLLCha =", cVLLCha ENDL
#endif


	cSLL1Neu = 0

	LOOP(Neu6, 1,4,1)
	LOOP(Neu5, 1,4,1)
	LOOP(All6, 1,6,1)
	LOOP(All5, 1,6,1)

        cSLL1Neu = cSLL1Neu + 
     &    2/(81.D0*Pi**2)*(GF**2*MW2**2*SW2*
     &        D0z(MASf2(All5,4),MASf2(All6,4),MNeu2(Neu5),
     &         MNeu2(Neu6))*MNeu(Neu5)*MNeu(Neu6)*UASf(All5,6,4)*
     &        UASf(All6,6,4)*UASfC(All5,2,4)*UASfC(All6,2,4)*
     &        ZNeu(Neu6,1)*
     &        (ZNeu(Neu5,1)**2*
     &           (2*SW2*ZNeu(Neu6,1) - 3*CW*SW*ZNeu(Neu6,2)) + 
     &          9*(CW2*ZNeu(Neu5,2)**2*ZNeu(Neu6,1) + 
     &             ZNeu(Neu5,1)*ZNeu(Neu5,2)*
     &              (-(CW*SW*ZNeu(Neu6,1)) + CW2*ZNeu(Neu6,2)))))/
     &      CW2**2

	ENDLOOP(All5)
	ENDLOOP(All6)
	ENDLOOP(Neu5)
	ENDLOOP(Neu6)

#ifdef DETAILED_DEBUG
	DCONST "cSLL1Neu =", cSLL1Neu ENDL
#endif

	cLR2Neu = 0

	LOOP(Neu6, 1,4,1)
	LOOP(Neu5, 1,4,1)
	LOOP(All6, 1,6,1)
	LOOP(All5, 1,6,1)

	dup1 = CW*SW*ZNeuC(Neu6,1) - 3*CW2*ZNeuC(Neu6,2)

        cLR2Neu = cLR2Neu + 
     &    8/(81.D0*Pi**2)*(GF**2*MW2**2*SW2*
     &        D00z(MASf2(All5,4),MASf2(All6,4),MNeu2(Neu5),
     &         MNeu2(Neu6))*
     &        (UASf(All5,6,4)*UASf(All6,3,4)*UASfC(All5,2,4)*
     &           UASfC(All6,5,4)*ZNeu(Neu6,1)*
     &           ((3*ZNeu(Neu5,2)*
     &                 (-2*CW*SW*ZNeuC(Neu5,1) + 
     &                   3*CW2*ZNeuC(Neu5,2)) + 
     &                ZNeu(Neu5,1)*
     &                 (2*SW2*ZNeuC(Neu5,1) - 
     &                   3*CW*SW*ZNeuC(Neu5,2)))*ZNeuC(Neu6,1) - 
     &             3*(CW*SW*ZNeu(Neu5,1) - 3*CW2*ZNeu(Neu5,2))*
     &              ZNeuC(Neu5,1)*ZNeuC(Neu6,2)) + 
     &          UASf(All5,3,4)*UASf(All6,6,4)*UASfC(All5,5,4)*
     &           UASfC(All6,2,4)*ZNeuC(Neu5,1)*
     &           (-3*dup1*ZNeu(Neu5,2)*ZNeu(Neu6,1) + 
     &             ZNeu(Neu5,1)*
     &              (-3*dup1*ZNeu(Neu6,2) + 
     &                2*ZNeu(Neu6,1)*
     &                 (SW2*ZNeuC(Neu6,1) - 3*CW*SW*ZNeuC(Neu6,2)))
     &             )))/CW2**2

	ENDLOOP(All5)
	ENDLOOP(All6)
	ENDLOOP(Neu5)
	ENDLOOP(Neu6)

#ifdef DETAILED_DEBUG
	DCONST "cLR2Neu =", cLR2Neu ENDL
#endif

	cSRR1Neu = 0

	LOOP(Neu6, 1,4,1)
	LOOP(Neu5, 1,4,1)
	LOOP(All6, 1,6,1)
	LOOP(All5, 1,6,1)

        cSRR1Neu = cSRR1Neu + 
     &    2/(81.D0*Pi**2)*(GF**2*MW2**2*SW2*
     &        D0z(MASf2(All5,4),MASf2(All6,4),MNeu2(Neu5),
     &         MNeu2(Neu6))*MNeu(Neu5)*MNeu(Neu6)*UASf(All5,3,4)*
     &        UASf(All6,3,4)*UASfC(All5,5,4)*UASfC(All6,5,4)*
     &        ZNeuC(Neu5,1)*
     &        ((2*SW2*ZNeuC(Neu5,1) - 3*CW*SW*ZNeuC(Neu5,2))*
     &           ZNeuC(Neu6,1)**2 + 
     &          9*((-(CW*SW*ZNeuC(Neu5,1)) + CW2*ZNeuC(Neu5,2))*
     &              ZNeuC(Neu6,1)*ZNeuC(Neu6,2) + 
     &             CW2*ZNeuC(Neu5,1)*ZNeuC(Neu6,2)**2)))/CW2**2

	ENDLOOP(All5)
	ENDLOOP(All6)
	ENDLOOP(Neu5)
	ENDLOOP(Neu6)

#ifdef DETAILED_DEBUG
	DCONST "cSRR1Neu =", cSRR1Neu ENDL
#endif

	cVLLNeu = 0

	LOOP(Neu6, 1,4,1)
	LOOP(Neu5, 1,4,1)
	LOOP(All6, 1,6,1)
	LOOP(All5, 1,6,1)

	dup2 = -(CW*SW*ZNeuC(Neu5,1)) + 3*CW2*ZNeuC(Neu5,2)

	dup3 = SW2*ZNeuC(Neu5,1) - 3*CW*SW*ZNeuC(Neu5,2)

	dup4 = CW*SW*ZNeuC(Neu6,1) - 3*CW2*ZNeuC(Neu6,2)

	dup5 = SW2*ZNeuC(Neu6,1) - 3*CW*SW*ZNeuC(Neu6,2)

        dup6 = SW2*ZNeuC(Neu6,1)**2 - 
     &    6*CW*SW*ZNeuC(Neu6,1)*ZNeuC(Neu6,2) + 
     &    9*CW2*ZNeuC(Neu6,2)**2

        cVLLNeu = cVLLNeu - 
     &    1/(162.D0*Pi**2)*(GF**2*MW2**2*UASf(All5,3,4)*
     &        UASf(All6,3,4)*UASfC(All5,2,4)*UASfC(All6,2,4)*
     &        (2*D00z(MASf2(All5,4),MASf2(All6,4),MNeu2(Neu5),
     &            MNeu2(Neu6))*
     &           (3*ZNeu(Neu5,2)*
     &              (3*CW2*ZNeu(Neu6,2)*
     &                 (dup5*ZNeuC(Neu5,1) - 3*dup4*ZNeuC(Neu5,2))+
     &                  ZNeu(Neu6,1)*
     &                 (dup2*SW2*ZNeuC(Neu6,1) + 
     &                   3*CW2*dup3*ZNeuC(Neu6,2))) + 
     &             ZNeu(Neu5,1)*
     &              (SW2*ZNeu(Neu6,1)*
     &                 (dup5*ZNeuC(Neu5,1) - 3*dup4*ZNeuC(Neu5,2))+
     &                  3*ZNeu(Neu6,2)*
     &                 (dup2*SW2*ZNeuC(Neu6,1) + 
     &                   3*CW2*dup3*ZNeuC(Neu6,2)))) + 
     &          D0z(MASf2(All5,4),MASf2(All6,4),MNeu2(Neu5),
     &            MNeu2(Neu6))*MNeu(Neu5)*MNeu(Neu6)*
     &           (dup6*(SW2*ZNeu(Neu5,1)**2 + 
     &                9*CW2*ZNeu(Neu5,2)**2) - 
     &             6*ZNeu(Neu5,1)*ZNeu(Neu5,2)*
     &              (-6*CW2*SW2*ZNeuC(Neu6,1)*ZNeuC(Neu6,2) + 
     &                CW*SW*
     &                 (SW2*ZNeuC(Neu6,1)**2 + 
     &                   9*CW2*ZNeuC(Neu6,2)**2)))))/CW2**2

	ENDLOOP(All5)
	ENDLOOP(All6)
	ENDLOOP(Neu5)
	ENDLOOP(Neu6)

#ifdef DETAILED_DEBUG
	DCONST "cVLLNeu =", cVLLNeu ENDL
#endif

	cLR1Neu = 0

	LOOP(Neu6, 1,4,1)
	LOOP(Neu5, 1,4,1)
	LOOP(All6, 1,6,1)
	LOOP(All5, 1,6,1)

	dup7 = CW*SW*ZNeuC(Neu6,1) - 3*CW2*ZNeuC(Neu6,2)

	dup8 = SW2*ZNeuC(Neu6,1) - 3*CW*SW*ZNeuC(Neu6,2)

        cLR1Neu = cLR1Neu + 
     &    2/(81.D0*Pi**2)*(GF**2*MW2**2*SW2*
     &        (2*D00z(MASf2(All5,4),MASf2(All6,4),MNeu2(Neu5),
     &            MNeu2(Neu6))*
     &           (UASf(All5,6,4)*UASf(All6,3,4)*UASfC(All5,5,4)*
     &              UASfC(All6,2,4) + 
     &             UASf(All5,3,4)*UASf(All6,6,4)*UASfC(All5,2,4)*
     &              UASfC(All6,5,4))*
     &           (dup8*ZNeu(Neu5,1) - 3*dup7*ZNeu(Neu5,2))*
     &           ZNeu(Neu6,1)*ZNeuC(Neu5,1) + 
     &          D0z(MASf2(All5,4),MASf2(All6,4),MNeu2(Neu5),
     &            MNeu2(Neu6))*MNeu(Neu5)*MNeu(Neu6)*
     &           (UASf(All5,6,4)*UASf(All6,3,4)*UASfC(All5,5,4)*
     &              UASfC(All6,2,4)*ZNeu(Neu6,1)*ZNeuC(Neu5,1)*
     &              (-3*ZNeu(Neu6,2)*
     &                 (CW*SW*ZNeuC(Neu5,1) - 3*CW2*ZNeuC(Neu5,2))+
     &                  ZNeu(Neu6,1)*
     &                 (SW2*ZNeuC(Neu5,1) - 3*CW*SW*ZNeuC(Neu5,2)))
     &               + UASf(All5,3,4)*UASf(All6,6,4)*
     &              UASfC(All5,2,4)*UASfC(All6,5,4)*ZNeu(Neu5,1)*
     &              (dup8*ZNeu(Neu5,1) - 3*dup7*ZNeu(Neu5,2))*
     &              ZNeuC(Neu6,1))))/CW2**2

	ENDLOOP(All5)
	ENDLOOP(All6)
	ENDLOOP(Neu5)
	ENDLOOP(Neu6)

#ifdef DETAILED_DEBUG
	DCONST "cLR1Neu =", cLR1Neu ENDL
#endif

	cVRRNeu = 0

	LOOP(Neu6, 1,4,1)
	LOOP(Neu5, 1,4,1)
	LOOP(All6, 1,6,1)
	LOOP(All5, 1,6,1)

        cVRRNeu = cVRRNeu - 
     &    8/(81.D0*Pi**2)*(GF**2*MW2**2*SW2**2*UASf(All5,6,4)*
     &        UASf(All6,6,4)*UASfC(All5,5,4)*UASfC(All6,5,4)*
     &        ZNeu(Neu6,1)*ZNeuC(Neu5,1)*
     &        (D0z(MASf2(All5,4),MASf2(All6,4),MNeu2(Neu5),
     &            MNeu2(Neu6))*MNeu(Neu5)*MNeu(Neu6)*ZNeu(Neu6,1)*
     &           ZNeuC(Neu5,1) + 
     &          2*D00z(MASf2(All5,4),MASf2(All6,4),MNeu2(Neu5),
     &            MNeu2(Neu6))*ZNeu(Neu5,1)*ZNeuC(Neu6,1)))/CW2**2

	ENDLOOP(All5)
	ENDLOOP(All6)
	ENDLOOP(Neu5)
	ENDLOOP(Neu6)

#ifdef DETAILED_DEBUG
	DCONST "cVRRNeu =", cVRRNeu ENDL
#endif

	cSRR2Neu = 0

	LOOP(Neu6, 1,4,1)
	LOOP(Neu5, 1,4,1)
	LOOP(All6, 1,6,1)
	LOOP(All5, 1,6,1)

        cSRR2Neu = cSRR2Neu + 
     &    1/(54.D0*Pi**2)*(GF**2*MW2**2*SW2*
     &        D0z(MASf2(All5,4),MASf2(All6,4),MNeu2(Neu5),
     &         MNeu2(Neu6))*MNeu(Neu5)*MNeu(Neu6)*UASf(All5,3,4)*
     &        UASf(All6,3,4)*UASfC(All5,5,4)*UASfC(All6,5,4)*
     &        ZNeuC(Neu5,1)*
     &        (CW*SW*ZNeuC(Neu6,1) - 3*CW2*ZNeuC(Neu6,2))*
     &        (ZNeuC(Neu5,2)*ZNeuC(Neu6,1) - 
     &          ZNeuC(Neu5,1)*ZNeuC(Neu6,2)))/CW2**2

	ENDLOOP(All5)
	ENDLOOP(All6)
	ENDLOOP(Neu5)
	ENDLOOP(Neu6)

#ifdef DETAILED_DEBUG
	DCONST "cSRR2Neu =", cSRR2Neu ENDL
#endif

	cSLL2Neu = 0

	LOOP(Neu6, 1,4,1)
	LOOP(Neu5, 1,4,1)
	LOOP(All6, 1,6,1)
	LOOP(All5, 1,6,1)

        cSLL2Neu = cSLL2Neu + 
     &    1/(54.D0*Pi**2)*(GF**2*MW2**2*SW2*
     &        D0z(MASf2(All5,4),MASf2(All6,4),MNeu2(Neu5),
     &         MNeu2(Neu6))*MNeu(Neu5)*MNeu(Neu6)*UASf(All5,6,4)*
     &        UASf(All6,6,4)*UASfC(All5,2,4)*UASfC(All6,2,4)*
     &        (CW*SW*ZNeu(Neu5,1) - 3*CW2*ZNeu(Neu5,2))*
     &        ZNeu(Neu6,1)*
     &        (-(ZNeu(Neu5,2)*ZNeu(Neu6,1)) + 
     &          ZNeu(Neu5,1)*ZNeu(Neu6,2)))/CW2**2

	ENDLOOP(All5)
	ENDLOOP(All6)
	ENDLOOP(Neu5)
	ENDLOOP(Neu6)

#ifdef DETAILED_DEBUG
	DCONST "cSLL2Neu =", cSLL2Neu ENDL
#endif


	cSLL1Glu = 0

	LOOP(All6, 1,6,1)
	LOOP(All5, 1,6,1)

        cSLL1Glu = cSLL1Glu - 
     &    37/9.D0*(M_3**2*asMT**2*
     &       D0z(MGl2,MGl2,MASf2(All5,4),MASf2(All6,4))*
     &       UASf(All5,6,4)*UASf(All6,6,4)*UASfC(All5,2,4)*
     &       UASfC(All6,2,4))

	ENDLOOP(All5)
	ENDLOOP(All6)

	LOOP(Neu5, 1,4,1)
	LOOP(All6, 1,6,1)
	LOOP(All5, 1,6,1)

        cSLL1Glu = cSLL1Glu + 
     &    2/(27.D0*Pi*sqrt2)*
     &     (M_3*asMT*GF*MW2*
     &        D0z(MGl2,MASf2(All5,4),MASf2(All6,4),MNeu2(Neu5))*
     &        MNeu(Neu5)*UASf(All5,6,4)*UASf(All6,6,4)*
     &        UASfC(All5,2,4)*UASfC(All6,2,4)*
     &        (9*CW2*SW*ZNeu(Neu5,2)**2 + 
     &          SW2*(19*SW*ZNeu(Neu5,1)**2 - 
     &             48*CW*ZNeu(Neu5,1)*ZNeu(Neu5,2))))/(CW2*SW)

	ENDLOOP(All5)
	ENDLOOP(All6)
	ENDLOOP(Neu5)

#ifdef DETAILED_DEBUG
	DCONST "cSLL1Glu =", cSLL1Glu ENDL
#endif

	cLR2Glu = 0

	LOOP(All6, 1,6,1)
	LOOP(All5, 1,6,1)

        cLR2Glu = cLR2Glu + 
     &    2/9.D0*(M_3*M_3C*asMT**2*
     &        (-21*MGl2*D0z(MGl2,MGl2,MASf2(All5,4),MASf2(All6,4))*
     &           (UASf(All5,6,4)*UASf(All6,3,4)*UASfC(All5,5,4)*
     &              UASfC(All6,2,4) + 
     &             UASf(All5,3,4)*UASf(All6,6,4)*UASfC(All5,2,4)*
     &              UASfC(All6,5,4)) + 
     &          2*D00z(MGl2,MGl2,MASf2(All5,4),MASf2(All6,4))*
     &           (UASf(All5,3,4)*UASf(All6,6,4)*
     &              (11*UASfC(All5,5,4)*UASfC(All6,2,4) + 
     &                6*UASfC(All5,2,4)*UASfC(All6,5,4)) + 
     &             UASf(All5,6,4)*UASf(All6,3,4)*
     &              (6*UASfC(All5,5,4)*UASfC(All6,2,4) + 
     &                11*UASfC(All5,2,4)*UASfC(All6,5,4)))))/MGl2

	ENDLOOP(All5)
	ENDLOOP(All6)

	LOOP(Neu5, 1,4,1)
	LOOP(All6, 1,6,1)
	LOOP(All5, 1,6,1)

	dup1 = -2*SW*ZNeuC(Neu5,1) + 3*CW*ZNeuC(Neu5,2)

	dup2 = SW*ZNeuC(Neu5,1) + 3*CW*ZNeuC(Neu5,2)

	dup3 = CW*SW2*ZNeuC(Neu5,1) + 3*CW2*SW*ZNeuC(Neu5,2)

        cLR2Glu = cLR2Glu - 
     &    4/(27.D0*Pi*sqrt2)*
     &     (asMT*GF*MW2*(D00z(MGl2,MASf2(All5,4),MASf2(All6,4),
     &            MNeu2(Neu5))*
     &           (UASf(All5,6,4)*UASf(All6,3,4)*
     &              (UASfC(All5,2,4)*UASfC(All6,5,4)*
     &                 (dup2*SW2*ZNeu(Neu5,1) + 
     &                   3*dup3*ZNeu(Neu5,2)) + 
     &                6*SW2*UASfC(All5,5,4)*UASfC(All6,2,4)*
     &                 (dup1*ZNeu(Neu5,1) + 
     &                   3*CW*ZNeu(Neu5,2)*ZNeuC(Neu5,1))) + 
     &             UASf(All5,3,4)*UASf(All6,6,4)*
     &              (UASfC(All5,5,4)*UASfC(All6,2,4)*
     &                 (dup2*SW2*ZNeu(Neu5,1) + 
     &                   3*dup3*ZNeu(Neu5,2)) + 
     &                6*SW2*UASfC(All5,2,4)*UASfC(All6,5,4)*
     &                 (dup1*ZNeu(Neu5,1) + 
     &                   3*CW*ZNeu(Neu5,2)*ZNeuC(Neu5,1)))) - 
     &          3*SW2*D0z(MGl2,MASf2(All5,4),MASf2(All6,4),
     &            MNeu2(Neu5))*MNeu(Neu5)*
     &           (UASf(All5,6,4)*UASf(All6,3,4)*UASfC(All5,5,4)*
     &              UASfC(All6,2,4) + 
     &             UASf(All5,3,4)*UASf(All6,6,4)*UASfC(All5,2,4)*
     &              UASfC(All6,5,4))*
     &           (M_3C*(SW*ZNeu(Neu5,1)**2 - 
     &                3*CW*ZNeu(Neu5,1)*ZNeu(Neu5,2)) + 
     &             M_3*ZNeuC(Neu5,1)*
     &              (SW*ZNeuC(Neu5,1) - 3*CW*ZNeuC(Neu5,2)))))/
     &      (CW2*SW)

	ENDLOOP(All5)
	ENDLOOP(All6)
	ENDLOOP(Neu5)

#ifdef DETAILED_DEBUG
	DCONST "cLR2Glu =", cLR2Glu ENDL
#endif

	cSRR1Glu = 0

	LOOP(All6, 1,6,1)
	LOOP(All5, 1,6,1)

        cSRR1Glu = cSRR1Glu - 
     &    37/9.D0*(M_3C**2*asMT**2*
     &       D0z(MGl2,MGl2,MASf2(All5,4),MASf2(All6,4))*
     &       UASf(All5,3,4)*UASf(All6,3,4)*UASfC(All5,5,4)*
     &       UASfC(All6,5,4))

	ENDLOOP(All5)
	ENDLOOP(All6)

	LOOP(Neu5, 1,4,1)
	LOOP(All6, 1,6,1)
	LOOP(All5, 1,6,1)

        cSRR1Glu = cSRR1Glu + 
     &    2/(27.D0*Pi*sqrt2)*
     &     (M_3C*asMT*GF*MW2*
     &        D0z(MGl2,MASf2(All5,4),MASf2(All6,4),MNeu2(Neu5))*
     &        MNeu(Neu5)*UASf(All5,3,4)*UASf(All6,3,4)*
     &        UASfC(All5,5,4)*UASfC(All6,5,4)*
     &        (9*CW2*SW*ZNeuC(Neu5,2)**2 + 
     &          SW2*(19*SW*ZNeuC(Neu5,1)**2 - 
     &             48*CW*ZNeuC(Neu5,1)*ZNeuC(Neu5,2))))/(CW2*SW)

	ENDLOOP(All5)
	ENDLOOP(All6)
	ENDLOOP(Neu5)

#ifdef DETAILED_DEBUG
	DCONST "cSRR1Glu =", cSRR1Glu ENDL
#endif

	cVLLGlu = 0

	LOOP(All6, 1,6,1)
	LOOP(All5, 1,6,1)

        cVLLGlu = cVLLGlu + 
     &    1/9.D0*(M_3*M_3C*asMT**2*
     &        (-44*D00z(MGl2,MGl2,MASf2(All5,4),MASf2(All6,4)) - 
     &          4*MGl2*D0z(MGl2,MGl2,MASf2(All5,4),MASf2(All6,4)))*
     &        UASf(All5,3,4)*UASf(All6,3,4)*UASfC(All5,2,4)*
     &        UASfC(All6,2,4))/MGl2

	ENDLOOP(All5)
	ENDLOOP(All6)

	LOOP(Neu5, 1,4,1)
	LOOP(All6, 1,6,1)
	LOOP(All5, 1,6,1)

        cVLLGlu = cVLLGlu - 
     &    2/(27.D0*Pi*sqrt2)*
     &     (asMT*GF*MW2*UASf(All5,3,4)*UASf(All6,3,4)*
     &        UASfC(All5,2,4)*UASfC(All6,2,4)*
     &        (4*D00z(MGl2,MASf2(All5,4),MASf2(All6,4),
     &            MNeu2(Neu5))*
     &           (SW2*ZNeu(Neu5,1)*
     &              (SW*ZNeuC(Neu5,1) - 3*CW*ZNeuC(Neu5,2)) - 
     &             3*ZNeu(Neu5,2)*
     &              (CW*SW2*ZNeuC(Neu5,1) - 3*CW2*SW*ZNeuC(Neu5,2))
     &             ) + D0z(MGl2,MASf2(All5,4),MASf2(All6,4),
     &            MNeu2(Neu5))*MNeu(Neu5)*
     &           (9*CW2*SW*
     &              (M_3C*ZNeu(Neu5,2)**2 + 
     &                M_3*ZNeuC(Neu5,2)**2) + 
     &             SW2*(SW*
     &                 (M_3C*ZNeu(Neu5,1)**2 + 
     &                   M_3*ZNeuC(Neu5,1)**2) - 
     &                6*CW*
     &                 (M_3C*ZNeu(Neu5,1)*ZNeu(Neu5,2) + 
     &                   M_3*ZNeuC(Neu5,1)*ZNeuC(Neu5,2))))))/
     &      (CW2*SW)

	ENDLOOP(All5)
	ENDLOOP(All6)
	ENDLOOP(Neu5)

#ifdef DETAILED_DEBUG
	DCONST "cVLLGlu =", cVLLGlu ENDL
#endif

	cLR1Glu = 0

	LOOP(All6, 1,6,1)
	LOOP(All5, 1,6,1)

        cLR1Glu = cLR1Glu + 
     &    1/9.D0*(M_3*M_3C*asMT**2*
     &        (MGl2*D0z(MGl2,MGl2,MASf2(All5,4),MASf2(All6,4))*
     &           (UASf(All5,6,4)*UASf(All6,3,4)*UASfC(All5,5,4)*
     &              UASfC(All6,2,4) + 
     &             UASf(All5,3,4)*UASf(All6,6,4)*UASfC(All5,2,4)*
     &              UASfC(All6,5,4)) + 
     &          10*D00z(MGl2,MGl2,MASf2(All5,4),MASf2(All6,4))*
     &           (UASf(All5,6,4)*UASf(All6,3,4)*
     &              (2*UASfC(All5,5,4)*UASfC(All6,2,4) - 
     &                3*UASfC(All5,2,4)*UASfC(All6,5,4)) + 
     &             UASf(All5,3,4)*UASf(All6,6,4)*
     &              (-3*UASfC(All5,5,4)*UASfC(All6,2,4) + 
     &                2*UASfC(All5,2,4)*UASfC(All6,5,4)))))/MGl2

	ENDLOOP(All5)
	ENDLOOP(All6)

	LOOP(Neu5, 1,4,1)
	LOOP(All6, 1,6,1)
	LOOP(All5, 1,6,1)

	dup4 = -2*SW*ZNeuC(Neu5,1) + 3*CW*ZNeuC(Neu5,2)

	dup5 = SW*ZNeuC(Neu5,1) + 3*CW*ZNeuC(Neu5,2)

	dup6 = CW*SW2*ZNeuC(Neu5,1) + 3*CW2*SW*ZNeuC(Neu5,2)

        cLR1Glu = cLR1Glu - 
     &    2/(27.D0*Pi*sqrt2)*
     &     (asMT*GF*MW2*(D00z(MGl2,MASf2(All5,4),MASf2(All6,4),
     &            MNeu2(Neu5))*
     &           (UASf(All5,6,4)*UASf(All6,3,4)*
     &              (3*UASfC(All5,2,4)*UASfC(All6,5,4)*
     &                 (dup5*SW2*ZNeu(Neu5,1) + 
     &                   3*dup6*ZNeu(Neu5,2)) + 
     &                2*SW2*UASfC(All5,5,4)*UASfC(All6,2,4)*
     &                 (dup4*ZNeu(Neu5,1) + 
     &                   3*CW*ZNeu(Neu5,2)*ZNeuC(Neu5,1))) + 
     &             UASf(All5,3,4)*UASf(All6,6,4)*
     &              (3*UASfC(All5,5,4)*UASfC(All6,2,4)*
     &                 (dup5*SW2*ZNeu(Neu5,1) + 
     &                   3*dup6*ZNeu(Neu5,2)) + 
     &                2*SW2*UASfC(All5,2,4)*UASfC(All6,5,4)*
     &                 (dup4*ZNeu(Neu5,1) + 
     &                   3*CW*ZNeu(Neu5,2)*ZNeuC(Neu5,1)))) - 
     &          SW2*D0z(MGl2,MASf2(All5,4),MASf2(All6,4),
     &            MNeu2(Neu5))*MNeu(Neu5)*
     &           (UASf(All5,6,4)*UASf(All6,3,4)*UASfC(All5,5,4)*
     &              UASfC(All6,2,4) + 
     &             UASf(All5,3,4)*UASf(All6,6,4)*UASfC(All5,2,4)*
     &              UASfC(All6,5,4))*
     &           (M_3C*(SW*ZNeu(Neu5,1)**2 - 
     &                3*CW*ZNeu(Neu5,1)*ZNeu(Neu5,2)) + 
     &             M_3*ZNeuC(Neu5,1)*
     &              (SW*ZNeuC(Neu5,1) - 3*CW*ZNeuC(Neu5,2)))))/
     &      (CW2*SW)

	ENDLOOP(All5)
	ENDLOOP(All6)
	ENDLOOP(Neu5)

#ifdef DETAILED_DEBUG
	DCONST "cLR1Glu =", cLR1Glu ENDL
#endif

	cVRRGlu = 0

	LOOP(All6, 1,6,1)
	LOOP(All5, 1,6,1)

        cVRRGlu = cVRRGlu + 
     &    1/9.D0*(M_3*M_3C*asMT**2*
     &        (-44*D00z(MGl2,MGl2,MASf2(All5,4),MASf2(All6,4)) - 
     &          4*MGl2*D0z(MGl2,MGl2,MASf2(All5,4),MASf2(All6,4)))*
     &        UASf(All5,6,4)*UASf(All6,6,4)*UASfC(All5,5,4)*
     &        UASfC(All6,5,4))/MGl2

	ENDLOOP(All5)
	ENDLOOP(All6)

	LOOP(Neu5, 1,4,1)
	LOOP(All6, 1,6,1)
	LOOP(All5, 1,6,1)

        cVRRGlu = cVRRGlu - 
     &    8/(27.D0*Pi*sqrt2)*
     &     (asMT*GF*MW2*SW2*UASf(All5,6,4)*UASf(All6,6,4)*
     &        UASfC(All5,5,4)*UASfC(All6,5,4)*
     &        (4*D00z(MGl2,MASf2(All5,4),MASf2(All6,4),
     &            MNeu2(Neu5))*ZNeu(Neu5,1)*ZNeuC(Neu5,1) + 
     &          D0z(MGl2,MASf2(All5,4),MASf2(All6,4),MNeu2(Neu5))*
     &           MNeu(Neu5)*
     &           (M_3C*ZNeu(Neu5,1)**2 + M_3*ZNeuC(Neu5,1)**2))
     &        )/CW2

	ENDLOOP(All5)
	ENDLOOP(All6)
	ENDLOOP(Neu5)

#ifdef DETAILED_DEBUG
	DCONST "cVRRGlu =", cVRRGlu ENDL
#endif

	cSRR2Glu = 0

	LOOP(All6, 1,6,1)
	LOOP(All5, 1,6,1)

        cSRR2Glu = cSRR2Glu + 
     &    1/12.D0*(M_3C**2*asMT**2*
     &       D0z(MGl2,MGl2,MASf2(All5,4),MASf2(All6,4))*
     &       UASf(All5,3,4)*UASf(All6,3,4)*UASfC(All5,5,4)*
     &       UASfC(All6,5,4))

	ENDLOOP(All5)
	ENDLOOP(All6)

	LOOP(Neu5, 1,4,1)
	LOOP(All6, 1,6,1)
	LOOP(All5, 1,6,1)

        cSRR2Glu = cSRR2Glu + 
     &    1/(18.D0*Pi*sqrt2)*
     &     (M_3C*asMT*GF*MW2*
     &        D0z(MGl2,MASf2(All5,4),MASf2(All6,4),MNeu2(Neu5))*
     &        MNeu(Neu5)*UASf(All5,3,4)*UASf(All6,3,4)*
     &        UASfC(All5,5,4)*UASfC(All6,5,4)*
     &        (SW2*ZNeuC(Neu5,1)**2 + 3*CW2*ZNeuC(Neu5,2)**2))/CW2

	ENDLOOP(All5)
	ENDLOOP(All6)
	ENDLOOP(Neu5)

#ifdef DETAILED_DEBUG
	DCONST "cSRR2Glu =", cSRR2Glu ENDL
#endif

	cSLL2Glu = 0

	LOOP(All6, 1,6,1)
	LOOP(All5, 1,6,1)

        cSLL2Glu = cSLL2Glu + 
     &    1/12.D0*(M_3**2*asMT**2*
     &       D0z(MGl2,MGl2,MASf2(All5,4),MASf2(All6,4))*
     &       UASf(All5,6,4)*UASf(All6,6,4)*UASfC(All5,2,4)*
     &       UASfC(All6,2,4))

	ENDLOOP(All5)
	ENDLOOP(All6)

	LOOP(Neu5, 1,4,1)
	LOOP(All6, 1,6,1)
	LOOP(All5, 1,6,1)

        cSLL2Glu = cSLL2Glu + 
     &    1/(18.D0*Pi*sqrt2)*
     &     (M_3*asMT*GF*MW2*
     &        D0z(MGl2,MASf2(All5,4),MASf2(All6,4),MNeu2(Neu5))*
     &        MNeu(Neu5)*UASf(All5,6,4)*UASf(All6,6,4)*
     &        UASfC(All5,2,4)*UASfC(All6,2,4)*
     &        (SW2*ZNeu(Neu5,1)**2 + 3*CW2*ZNeu(Neu5,2)**2))/CW2

	ENDLOOP(All5)
	ENDLOOP(All6)
	ENDLOOP(Neu5)

#ifdef DETAILED_DEBUG
	DCONST "cSLL2Glu =", cSLL2Glu ENDL
#endif

