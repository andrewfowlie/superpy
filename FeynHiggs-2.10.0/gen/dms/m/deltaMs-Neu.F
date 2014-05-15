
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

