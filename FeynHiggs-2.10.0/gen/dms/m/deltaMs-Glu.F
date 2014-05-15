
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

