
	CSLNeu = 0

	LOOP(Sfe5, 1,2,1)
	LOOP(Neu6, 1,4,1)
	LOOP(Neu5, 1,4,1)
	LOOP(All5, 1,6,1)

	dup1 = 2*SW*ZNeu(Neu6,1) + CW*ZNeu(Neu6,2)

	dup2 = 2*CW*SW2*ZNeu(Neu6,1) + CW2*SW*ZNeu(Neu6,2)

	dup3 = 2*SW*ZNeuC(Neu6,1) + CW*ZNeuC(Neu6,2)

	dup4 = 2*CW*SW2*ZNeuC(Neu6,1) + CW2*SW*ZNeuC(Neu6,2)

        CSLNeu = CSLNeu + 
     &    1/18.D0*(SW2*(4*D00z(MASf2(All5,bTR),MNeu2(Neu5),
     &            MNeu2(Neu6),MSf2(Sfe5,2,2))*USf(Sfe5,2,2,2)*
     &           USfC(Sfe5,1,2,2)*
     &           (2*UASfC(All5,5,bTR)*
     &              (CB2*MW2*UASf(All5,3,bTR)*
     &                 (ZNeu(Neu5,1)*
     &                    (dup1*SW2*ZNeuC(Neu5,1) - 
     &                      3*dup2*ZNeuC(Neu5,2)) + 
     &                   ZNeu(Neu5,2)*ZNeu(Neu6,1)*
     &                    (CW*SW2*ZNeuC(Neu5,1) - 
     &                      3*CW2*SW*ZNeuC(Neu5,2))) + 
     &                3*CB*MW*Mf(bTR,3)*UASf(All5,6,bTR)*
     &                 (dup2*ZNeu(Neu5,1) + 
     &                   CW2*SW*ZNeu(Neu5,2)*ZNeu(Neu6,1))*
     &                 ZNeuC(Neu5,3))*ZNeuC(Neu6,1) + 
     &             3*Mf(bTR,2)*UASfC(All5,2,bTR)*
     &              (CB*MW*UASf(All5,3,bTR)*
     &                 (CW2*ZNeu(Neu5,2)*ZNeu(Neu6,1)*
     &                    (SW*ZNeuC(Neu5,1) - 3*CW*ZNeuC(Neu5,2))+
     &                     ZNeu(Neu5,1)*
     &                    (dup2*ZNeuC(Neu5,1) - 
     &                      3*CW2*dup1*ZNeuC(Neu5,2))) + 
     &                3*CW2*Mf(bTR,3)*UASf(All5,6,bTR)*
     &                 (dup1*ZNeu(Neu5,1) + 
     &                   CW*ZNeu(Neu5,2)*ZNeu(Neu6,1))*
     &                 ZNeuC(Neu5,3))*ZNeuC(Neu6,3)) + 
     &          D0z(MASf2(All5,bTR),MNeu2(Neu5),MNeu2(Neu6),
     &            MSf2(Sfe5,2,2))*MNeu(Neu5)*
     &           (UASf(All5,6,bTR)*
     &              (2*UASfC(All5,5,bTR)*ZNeuC(Neu6,1)*
     &                 (2*CB2*MB*MW2*SW2*USf(Sfe5,2,2,2)*
     &                    USfC(Sfe5,1,2,2)*ZNeu(Neu5,1)*
     &                    (-3*CW*ZNeu(Neu5,2)*ZNeu(Neu6,1) + 
     &                      ZNeu(Neu5,1)*
     &                       (-2*SW*ZNeu(Neu6,1) + CW*ZNeu(Neu6,2))
     &                      ) + 
     &                   3*CB*MW*Mf(bTR,3)*MNeu(Neu6)*
     &                    USf(Sfe5,1,2,2)*USfC(Sfe5,2,2,2)*
     &                    ZNeuC(Neu5,3)*
     &                    (dup4*ZNeuC(Neu5,1) + 
     &                      CW2*SW*ZNeuC(Neu5,2)*ZNeuC(Neu6,1))) + 
     &                3*Mf(bTR,2)*UASfC(All5,2,bTR)*
     &                 (2*CB*MB*MW*USf(Sfe5,2,2,2)*
     &                    USfC(Sfe5,1,2,2)*ZNeu(Neu5,1)*
     &                    (-3*CW2*SW*ZNeu(Neu5,2)*ZNeu(Neu6,1) + 
     &                      ZNeu(Neu5,1)*
     &                       (-2*CW*SW2*ZNeu(Neu6,1) + 
     &                       CW2*SW*ZNeu(Neu6,2))) + 
     &                   3*CW2*Mf(bTR,3)*MNeu(Neu6)*
     &                    USf(Sfe5,1,2,2)*USfC(Sfe5,2,2,2)*
     &                    ZNeuC(Neu5,3)*
     &                    (dup3*ZNeuC(Neu5,1) + 
     &                      CW*ZNeuC(Neu5,2)*ZNeuC(Neu6,1)))*
     &                 ZNeuC(Neu6,3)) - 
     &             UASf(All5,3,bTR)*
     &              (3*MB*Mf(bTR,3)*USf(Sfe5,2,2,2)*
     &                 USfC(Sfe5,1,2,2)*ZNeu(Neu5,3)*
     &                 (2*CB*MW*UASfC(All5,5,bTR)*
     &                    (3*CW2*SW*ZNeu(Neu5,2)*ZNeu(Neu6,1) + 
     &                      ZNeu(Neu5,1)*
     &                       (2*CW*SW2*ZNeu(Neu6,1) - 
     &                       CW2*SW*ZNeu(Neu6,2)))*ZNeuC(Neu6,1) + 
     &                   3*CW2*Mf(bTR,2)*UASfC(All5,2,bTR)*
     &                    (3*CW*ZNeu(Neu5,2)*ZNeu(Neu6,1) + 
     &                      ZNeu(Neu5,1)*
     &                       (2*SW*ZNeu(Neu6,1) - CW*ZNeu(Neu6,2)))
     &                     *ZNeuC(Neu6,3)) + 
     &                MNeu(Neu6)*USf(Sfe5,1,2,2)*USfC(Sfe5,2,2,2)*
     &                 (2*CB2*MW2*UASfC(All5,5,bTR)*ZNeuC(Neu6,1)*
     &                    (-(dup3*SW2*ZNeuC(Neu5,1)**2) + 
     &                      3*CW2*SW*ZNeuC(Neu5,2)**2*
     &                       ZNeuC(Neu6,1) + 
     &                      ZNeuC(Neu5,1)*ZNeuC(Neu5,2)*
     &                       (5*CW*SW2*ZNeuC(Neu6,1) + 
     &                       3*CW2*SW*ZNeuC(Neu6,2))) + 
     &                   3*CB*MW*Mf(bTR,2)*UASfC(All5,2,bTR)*
     &                    (-(dup4*ZNeuC(Neu5,1)**2) + 
     &                      CW2*
     &                       (3*CW*ZNeuC(Neu5,2)**2*
     &                      ZNeuC(Neu6,1) + 
     &                       ZNeuC(Neu5,1)*ZNeuC(Neu5,2)*
     &                       (5*SW*ZNeuC(Neu6,1) + 
     &                       3*CW*ZNeuC(Neu6,2))))*ZNeuC(Neu6,3))))
     &          ))/(CB2*CW2**2*SW*CKM(3,3)*CKMC(3,2))

	ENDLOOP(All5)
	ENDLOOP(Neu5)
	ENDLOOP(Neu6)
	ENDLOOP(Sfe5)

#ifdef DETAILED_DEBUG
	DCONST "CSLNeu =", CSLNeu ENDL
#endif

	CPLNeu = 0

	LOOP(Sfe5, 1,2,1)
	LOOP(Neu6, 1,4,1)
	LOOP(Neu5, 1,4,1)
	LOOP(All5, 1,6,1)

	dup5 = 2*SW*ZNeu(Neu6,1) + CW*ZNeu(Neu6,2)

	dup6 = 2*CW*SW2*ZNeu(Neu6,1) + CW2*SW*ZNeu(Neu6,2)

	dup7 = 2*SW*ZNeuC(Neu6,1) + CW*ZNeuC(Neu6,2)

	dup8 = 2*CW*SW2*ZNeuC(Neu6,1) + CW2*SW*ZNeuC(Neu6,2)

        CPLNeu = CPLNeu - 
     &    1/18.D0*(SW2*(-4*D00z(MASf2(All5,bTR),MNeu2(Neu5),
     &            MNeu2(Neu6),MSf2(Sfe5,2,2))*USf(Sfe5,2,2,2)*
     &           USfC(Sfe5,1,2,2)*
     &           (2*UASfC(All5,5,bTR)*
     &              (CB2*MW2*UASf(All5,3,bTR)*
     &                 (ZNeu(Neu5,1)*
     &                    (dup5*SW2*ZNeuC(Neu5,1) - 
     &                      3*dup6*ZNeuC(Neu5,2)) + 
     &                   ZNeu(Neu5,2)*ZNeu(Neu6,1)*
     &                    (CW*SW2*ZNeuC(Neu5,1) - 
     &                      3*CW2*SW*ZNeuC(Neu5,2))) + 
     &                3*CB*MW*Mf(bTR,3)*UASf(All5,6,bTR)*
     &                 (dup6*ZNeu(Neu5,1) + 
     &                   CW2*SW*ZNeu(Neu5,2)*ZNeu(Neu6,1))*
     &                 ZNeuC(Neu5,3))*ZNeuC(Neu6,1) + 
     &             3*Mf(bTR,2)*UASfC(All5,2,bTR)*
     &              (CB*MW*UASf(All5,3,bTR)*
     &                 (CW2*ZNeu(Neu5,2)*ZNeu(Neu6,1)*
     &                    (SW*ZNeuC(Neu5,1) - 3*CW*ZNeuC(Neu5,2))+
     &                     ZNeu(Neu5,1)*
     &                    (dup6*ZNeuC(Neu5,1) - 
     &                      3*CW2*dup5*ZNeuC(Neu5,2))) + 
     &                3*CW2*Mf(bTR,3)*UASf(All5,6,bTR)*
     &                 (dup5*ZNeu(Neu5,1) + 
     &                   CW*ZNeu(Neu5,2)*ZNeu(Neu6,1))*
     &                 ZNeuC(Neu5,3))*ZNeuC(Neu6,3)) + 
     &          D0z(MASf2(All5,bTR),MNeu2(Neu5),MNeu2(Neu6),
     &            MSf2(Sfe5,2,2))*MNeu(Neu5)*
     &           (UASf(All5,6,bTR)*
     &              (2*UASfC(All5,5,bTR)*ZNeuC(Neu6,1)*
     &                 (2*CB2*MB*MW2*SW2*USf(Sfe5,2,2,2)*
     &                    USfC(Sfe5,1,2,2)*ZNeu(Neu5,1)*
     &                    (dup5*ZNeu(Neu5,1) + 
     &                      CW*ZNeu(Neu5,2)*ZNeu(Neu6,1)) + 
     &                   3*CB*MW*Mf(bTR,3)*MNeu(Neu6)*
     &                    USf(Sfe5,1,2,2)*USfC(Sfe5,2,2,2)*
     &                    ZNeuC(Neu5,3)*
     &                    (dup8*ZNeuC(Neu5,1) + 
     &                      CW2*SW*ZNeuC(Neu5,2)*ZNeuC(Neu6,1))) + 
     &                3*Mf(bTR,2)*UASfC(All5,2,bTR)*
     &                 (2*CB*MB*MW*USf(Sfe5,2,2,2)*
     &                    USfC(Sfe5,1,2,2)*ZNeu(Neu5,1)*
     &                    (dup6*ZNeu(Neu5,1) + 
     &                      CW2*SW*ZNeu(Neu5,2)*ZNeu(Neu6,1)) + 
     &                   3*CW2*Mf(bTR,3)*MNeu(Neu6)*
     &                    USf(Sfe5,1,2,2)*USfC(Sfe5,2,2,2)*
     &                    ZNeuC(Neu5,3)*
     &                    (dup7*ZNeuC(Neu5,1) + 
     &                      CW*ZNeuC(Neu5,2)*ZNeuC(Neu6,1)))*
     &                 ZNeuC(Neu6,3)) + 
     &             UASf(All5,3,bTR)*
     &              (3*MB*Mf(bTR,3)*USf(Sfe5,2,2,2)*
     &                 USfC(Sfe5,1,2,2)*ZNeu(Neu5,3)*
     &                 (2*CB*MW*UASfC(All5,5,bTR)*
     &                    (dup6*ZNeu(Neu5,1) + 
     &                      CW2*SW*ZNeu(Neu5,2)*ZNeu(Neu6,1))*
     &                    ZNeuC(Neu6,1) + 
     &                   3*CW2*Mf(bTR,2)*UASfC(All5,2,bTR)*
     &                    (dup5*ZNeu(Neu5,1) + 
     &                      CW*ZNeu(Neu5,2)*ZNeu(Neu6,1))*
     &                    ZNeuC(Neu6,3)) + 
     &                MNeu(Neu6)*USf(Sfe5,1,2,2)*USfC(Sfe5,2,2,2)*
     &                 (2*CB2*MW2*UASfC(All5,5,bTR)*ZNeuC(Neu6,1)*
     &                    (dup7*SW2*ZNeuC(Neu5,1)**2 - 
     &                      3*CW2*SW*ZNeuC(Neu5,2)**2*
     &                       ZNeuC(Neu6,1) - 
     &                      ZNeuC(Neu5,1)*ZNeuC(Neu5,2)*
     &                       (5*CW*SW2*ZNeuC(Neu6,1) + 
     &                       3*CW2*SW*ZNeuC(Neu6,2))) + 
     &                   3*CB*MW*Mf(bTR,2)*UASfC(All5,2,bTR)*
     &                    (dup8*ZNeuC(Neu5,1)**2 + 
     &                      CW2*
     &                       (-3*CW*ZNeuC(Neu5,2)**2*
     &                      ZNeuC(Neu6,1) - 
     &                       ZNeuC(Neu5,1)*ZNeuC(Neu5,2)*
     &                       (5*SW*ZNeuC(Neu6,1) + 
     &                       3*CW*ZNeuC(Neu6,2))))*ZNeuC(Neu6,3))))
     &          ))/(CB2*CW2**2*SW*CKM(3,3)*CKMC(3,2))

	ENDLOOP(All5)
	ENDLOOP(Neu5)
	ENDLOOP(Neu6)
	ENDLOOP(Sfe5)

#ifdef DETAILED_DEBUG
	DCONST "CPLNeu =", CPLNeu ENDL
#endif

	CSRNeu = 0

	LOOP(Sfe5, 1,2,1)
	LOOP(Neu6, 1,4,1)
	LOOP(Neu5, 1,4,1)
	LOOP(All5, 1,6,1)

	dup9 = SW*ZNeu(Neu6,1) - 3*CW*ZNeu(Neu6,2)

	dup10 = 2*SW*ZNeuC(Neu6,1) + CW*ZNeuC(Neu6,2)

	dup11 = 2*CW*SW2*ZNeuC(Neu6,1) + CW2*SW*ZNeuC(Neu6,2)

        tmp1 = 4*D00z(MASf2(All5,bTR),MNeu2(Neu5),MNeu2(Neu6),
     &     MSf2(Sfe5,2,2))*USf(Sfe5,1,2,2)*USfC(Sfe5,2,2,2)*
     &    (2*UASf(All5,6,bTR)*ZNeu(Neu5,1)*
     &       (3*CB*MW*Mf(bTR,2)*UASfC(All5,5,bTR)*ZNeu(Neu6,3)*
     &          (dup11*ZNeuC(Neu5,1) + 
     &            CW2*SW*ZNeuC(Neu5,2)*ZNeuC(Neu6,1)) + 
     &         CB2*MW2*UASfC(All5,2,bTR)*
     &          (SW2*ZNeu(Neu6,1)*
     &             (dup10*ZNeuC(Neu5,1) + 
     &               CW*ZNeuC(Neu5,2)*ZNeuC(Neu6,1)) - 
     &            3*ZNeu(Neu6,2)*
     &             (dup11*ZNeuC(Neu5,1) + 
     &               CW2*SW*ZNeuC(Neu5,2)*ZNeuC(Neu6,1)))) + 
     &      3*Mf(bTR,3)*UASf(All5,3,bTR)*ZNeu(Neu5,3)*
     &       (3*CW2*Mf(bTR,2)*UASfC(All5,5,bTR)*ZNeu(Neu6,3)*
     &          (dup10*ZNeuC(Neu5,1) + 
     &            CW*ZNeuC(Neu5,2)*ZNeuC(Neu6,1)) + 
     &         CB*MW*UASfC(All5,2,bTR)*
     &          ((-3*CW2*ZNeu(Neu6,2)*
     &                (2*SW*ZNeuC(Neu5,1) + CW*ZNeuC(Neu5,2)) + 
     &               ZNeu(Neu6,1)*
     &                (2*CW*SW2*ZNeuC(Neu5,1) + 
     &                  CW2*SW*ZNeuC(Neu5,2)))*ZNeuC(Neu6,1) + 
     &            CW2*dup9*ZNeuC(Neu5,1)*ZNeuC(Neu6,2))))

        tmp1 = tmp1 + D0z(MASf2(All5,bTR),MNeu2(Neu5),
     &      MNeu2(Neu6),MSf2(Sfe5,2,2))*MNeu(Neu5)*
     &     (MNeu(Neu6)*USf(Sfe5,2,2,2)*USfC(Sfe5,1,2,2)*
     &        (3*Mf(bTR,3)*UASf(All5,3,bTR)*ZNeu(Neu5,3)*
     &           (CB*MW*UASfC(All5,2,bTR)*
     &              (CW2*dup9*ZNeu(Neu5,2)*ZNeu(Neu6,1) + 
     &                ZNeu(Neu5,1)*
     &                 (-5*CW2*SW*ZNeu(Neu6,1)*ZNeu(Neu6,2) + 
     &                   CW*
     &                    (2*SW2*ZNeu(Neu6,1)**2 - 
     &                      3*CW2*ZNeu(Neu6,2)**2))) + 
     &             3*CW2*Mf(bTR,2)*UASfC(All5,5,bTR)*
     &              (CW*ZNeu(Neu5,2)*ZNeu(Neu6,1) + 
     &                ZNeu(Neu5,1)*
     &                 (2*SW*ZNeu(Neu6,1) + CW*ZNeu(Neu6,2)))*
     &              ZNeu(Neu6,3)) + 
     &          2*UASf(All5,6,bTR)*ZNeu(Neu5,1)*
     &           (CB2*MW2*UASfC(All5,2,bTR)*
     &              (SW2*(2*SW*ZNeu(Neu5,1) + CW*ZNeu(Neu5,2))*
     &                 ZNeu(Neu6,1)**2 - 
     &                (5*CW*SW2*ZNeu(Neu5,1) + 
     &                   3*CW2*SW*ZNeu(Neu5,2))*ZNeu(Neu6,1)*
     &                 ZNeu(Neu6,2) - 
     &                3*CW2*SW*ZNeu(Neu5,1)*ZNeu(Neu6,2)**2) + 
     &             3*CB*MW*Mf(bTR,2)*UASfC(All5,5,bTR)*
     &              (CW2*SW*ZNeu(Neu5,2)*ZNeu(Neu6,1) + 
     &                ZNeu(Neu5,1)*
     &                 (2*CW*SW2*ZNeu(Neu6,1) + 
     &                   CW2*SW*ZNeu(Neu6,2)))*ZNeu(Neu6,3))) - 
     &       MB*USf(Sfe5,1,2,2)*USfC(Sfe5,2,2,2)*
     &        (3*Mf(bTR,3)*UASf(All5,6,bTR)*ZNeuC(Neu5,3)*
     &           (CB*MW*UASfC(All5,2,bTR)*
     &              ((3*CW2*ZNeu(Neu6,2)*
     &                    (-2*SW*ZNeuC(Neu5,1) + CW*ZNeuC(Neu5,2))+
     &                     ZNeu(Neu6,1)*
     &                    (2*CW*SW2*ZNeuC(Neu5,1) - 
     &                      CW2*SW*ZNeuC(Neu5,2)))*ZNeuC(Neu6,1) + 
     &                3*CW2*dup9*ZNeuC(Neu5,1)*ZNeuC(Neu6,2)) + 
     &             3*CW2*Mf(bTR,2)*UASfC(All5,5,bTR)*ZNeu(Neu6,3)*
     &              (-(CW*ZNeuC(Neu5,2)*ZNeuC(Neu6,1)) + 
     &                ZNeuC(Neu5,1)*
     &                 (2*SW*ZNeuC(Neu6,1) + 3*CW*ZNeuC(Neu6,2))))+
     &            UASf(All5,3,bTR)*
     &           (CB2*MW2*UASfC(All5,2,bTR)*
     &              ((-7*(CW*SW2*ZNeu(Neu6,1) - 
     &                      3*CW2*SW*ZNeu(Neu6,2))*ZNeuC(Neu5,1)*
     &                    ZNeuC(Neu5,2) + 
     &                   dup9*
     &                    (2*SW2*ZNeuC(Neu5,1)**2 + 
     &                      3*CW2*ZNeuC(Neu5,2)**2))*ZNeuC(Neu6,1)+
     &                  3*ZNeuC(Neu5,1)*
     &                 (-3*CW2*ZNeu(Neu6,2)*
     &                    (SW*ZNeuC(Neu5,1) - 3*CW*ZNeuC(Neu5,2))+
     &                     ZNeu(Neu6,1)*
     &                    (CW*SW2*ZNeuC(Neu5,1) - 
     &                      3*CW2*SW*ZNeuC(Neu5,2)))*ZNeuC(Neu6,2))
     &               + 3*CB*MW*Mf(bTR,2)*UASfC(All5,5,bTR)*
     &              ZNeu(Neu6,3)*
     &              (ZNeuC(Neu5,1)**2*
     &                 (2*CW*SW2*ZNeuC(Neu6,1) + 
     &                   3*CW2*SW*ZNeuC(Neu6,2)) + 
     &                CW2*(3*CW*ZNeuC(Neu5,2)**2*ZNeuC(Neu6,1) - 
     &                   ZNeuC(Neu5,1)*ZNeuC(Neu5,2)*
     &                    (7*SW*ZNeuC(Neu6,1) + 9*CW*ZNeuC(Neu6,2))
     &                   )))))

        CSRNeu = CSRNeu + 
     &    1/18.D0*(SW2*tmp1)/(CB2*CW2**2*SW*CKM(3,3)*CKMC(3,2))

	ENDLOOP(All5)
	ENDLOOP(Neu5)
	ENDLOOP(Neu6)
	ENDLOOP(Sfe5)

#ifdef DETAILED_DEBUG
	DCONST "CSRNeu =", CSRNeu ENDL
#endif

	CPRNeu = 0

	LOOP(Sfe5, 1,2,1)
	LOOP(Neu6, 1,4,1)
	LOOP(Neu5, 1,4,1)
	LOOP(All5, 1,6,1)

	dup12 = SW*ZNeu(Neu6,1) - 3*CW*ZNeu(Neu6,2)

	dup13 = 2*SW*ZNeuC(Neu5,1) + CW*ZNeuC(Neu5,2)

	dup14 = 2*CW*SW2*ZNeuC(Neu5,1) + CW2*SW*ZNeuC(Neu5,2)

	dup15 = 2*SW*ZNeuC(Neu6,1) + CW*ZNeuC(Neu6,2)

	dup16 = 2*CW*SW2*ZNeuC(Neu6,1) + CW2*SW*ZNeuC(Neu6,2)

        dup17 = (dup14*ZNeu(Neu6,1) - 
     &       3*CW2*dup13*ZNeu(Neu6,2))*ZNeuC(Neu6,1) + 
     &    CW2*dup12*ZNeuC(Neu5,1)*ZNeuC(Neu6,2)

        CPRNeu = CPRNeu + 
     &    1/18.D0*(SW2*(-4*D00z(MASf2(All5,bTR),MNeu2(Neu5),
     &            MNeu2(Neu6),MSf2(Sfe5,2,2))*USf(Sfe5,1,2,2)*
     &           USfC(Sfe5,2,2,2)*
     &           (3*Mf(bTR,3)*UASf(All5,3,bTR)*ZNeu(Neu5,3)*
     &              (CB*dup17*MW*UASfC(All5,2,bTR) + 
     &                3*CW2*Mf(bTR,2)*UASfC(All5,5,bTR)*
     &                 ZNeu(Neu6,3)*
     &                 (dup15*ZNeuC(Neu5,1) + 
     &                   CW*ZNeuC(Neu5,2)*ZNeuC(Neu6,1))) + 
     &             2*UASf(All5,6,bTR)*ZNeu(Neu5,1)*
     &              (3*CB*MW*Mf(bTR,2)*UASfC(All5,5,bTR)*
     &                 ZNeu(Neu6,3)*
     &                 (dup16*ZNeuC(Neu5,1) + 
     &                   CW2*SW*ZNeuC(Neu5,2)*ZNeuC(Neu6,1)) + 
     &                CB2*MW2*UASfC(All5,2,bTR)*
     &                 (SW2*ZNeu(Neu6,1)*
     &                    (dup15*ZNeuC(Neu5,1) + 
     &                      CW*ZNeuC(Neu5,2)*ZNeuC(Neu6,1)) - 
     &                   3*ZNeu(Neu6,2)*
     &                    (dup16*ZNeuC(Neu5,1) + 
     &                      CW2*SW*ZNeuC(Neu5,2)*ZNeuC(Neu6,1)))))+
     &            D0z(MASf2(All5,bTR),MNeu2(Neu5),MNeu2(Neu6),
     &            MSf2(Sfe5,2,2))*MNeu(Neu5)*
     &           (MNeu(Neu6)*USf(Sfe5,2,2,2)*USfC(Sfe5,1,2,2)*
     &              (3*Mf(bTR,3)*UASf(All5,3,bTR)*ZNeu(Neu5,3)*
     &                 (CB*MW*UASfC(All5,2,bTR)*
     &                    (CW2*dup12*ZNeu(Neu5,2)*ZNeu(Neu6,1) + 
     &                      ZNeu(Neu5,1)*
     &                       (-5*CW2*SW*ZNeu(Neu6,1)*
     &                      ZNeu(Neu6,2) + 
     &                       CW*
     &                       (2*SW2*ZNeu(Neu6,1)**2 - 
     &                       3*CW2*ZNeu(Neu6,2)**2))) + 
     &                   3*CW2*Mf(bTR,2)*UASfC(All5,5,bTR)*
     &                    (CW*ZNeu(Neu5,2)*ZNeu(Neu6,1) + 
     &                      ZNeu(Neu5,1)*
     &                       (2*SW*ZNeu(Neu6,1) + CW*ZNeu(Neu6,2)))
     &                     *ZNeu(Neu6,3)) + 
     &                2*UASf(All5,6,bTR)*ZNeu(Neu5,1)*
     &                 (CB2*MW2*UASfC(All5,2,bTR)*
     &                    (SW2*
     &                       (2*SW*ZNeu(Neu5,1) + CW*ZNeu(Neu5,2))*
     &                       ZNeu(Neu6,1)**2 - 
     &                      (5*CW*SW2*ZNeu(Neu5,1) + 
     &                       3*CW2*SW*ZNeu(Neu5,2))*ZNeu(Neu6,1)*
     &                       ZNeu(Neu6,2) - 
     &                      3*CW2*SW*ZNeu(Neu5,1)*ZNeu(Neu6,2)**2)+
     &                     3*CB*MW*Mf(bTR,2)*UASfC(All5,5,bTR)*
     &                    (CW2*SW*ZNeu(Neu5,2)*ZNeu(Neu6,1) + 
     &                      ZNeu(Neu5,1)*
     &                       (2*CW*SW2*ZNeu(Neu6,1) + 
     &                       CW2*SW*ZNeu(Neu6,2)))*ZNeu(Neu6,3)))+
     &               MB*USf(Sfe5,1,2,2)*USfC(Sfe5,2,2,2)*
     &              (3*Mf(bTR,3)*UASf(All5,6,bTR)*ZNeuC(Neu5,3)*
     &                 (CB*dup17*MW*UASfC(All5,2,bTR) + 
     &                   3*CW2*Mf(bTR,2)*UASfC(All5,5,bTR)*
     &                    ZNeu(Neu6,3)*
     &                    (dup15*ZNeuC(Neu5,1) + 
     &                      CW*ZNeuC(Neu5,2)*ZNeuC(Neu6,1))) + 
     &                UASf(All5,3,bTR)*
     &                 (CB2*MW2*UASfC(All5,2,bTR)*
     &                    ((-5*
     &                       (CW*SW2*ZNeu(Neu6,1) - 
     &                       3*CW2*SW*ZNeu(Neu6,2))*ZNeuC(Neu5,1)*
     &                       ZNeuC(Neu5,2) + 
     &                       dup12*
     &                       (2*SW2*ZNeuC(Neu5,1)**2 - 
     &                       3*CW2*ZNeuC(Neu5,2)**2))*ZNeuC(Neu6,1)
     &                        + ZNeuC(Neu5,1)*
     &                       (-3*CW2*ZNeu(Neu6,2)*
     &                       (SW*ZNeuC(Neu5,1) - 
     &                       3*CW*ZNeuC(Neu5,2)) + 
     &                       ZNeu(Neu6,1)*
     &                       (CW*SW2*ZNeuC(Neu5,1) - 
     &                       3*CW2*SW*ZNeuC(Neu5,2)))*ZNeuC(Neu6,2)
     &                      ) + 
     &                   3*CB*MW*Mf(bTR,2)*UASfC(All5,5,bTR)*
     &                    ZNeu(Neu6,3)*
     &                    (dup16*ZNeuC(Neu5,1)**2 + 
     &                      CW2*
     &                       (-3*CW*ZNeuC(Neu5,2)**2*
     &                      ZNeuC(Neu6,1) - 
     &                       ZNeuC(Neu5,1)*ZNeuC(Neu5,2)*
     &                       (5*SW*ZNeuC(Neu6,1) + 
     &                       3*CW*ZNeuC(Neu6,2)))))))))/
     &      (CB2*CW2**2*SW*CKM(3,3)*CKMC(3,2))

	ENDLOOP(All5)
	ENDLOOP(Neu5)
	ENDLOOP(Neu6)
	ENDLOOP(Sfe5)

#ifdef DETAILED_DEBUG
	DCONST "CPRNeu =", CPRNeu ENDL
#endif

	CALNeu = 0

	LOOP(Neu5, 1,4,1)
	LOOP(All5, 1,6,1)

	dup18 = CW*SW*ZNeuC(Neu5,1) - 3*CW2*ZNeuC(Neu5,2)

        CALNeu = CALNeu - 
     &    1/864.D0*((1 + 2*CW2)*
     &        (-MASf2(All5,bTR)**2 + 
     &          2*A0(MASf2(All5,bTR))*
     &           (MASf2(All5,bTR) - 2*MNeu2(Neu5)) + 
     &          MNeu2(Neu5)*(2*A0(MNeu2(Neu5)) + MNeu2(Neu5)))*
     &        (UASf(All5,3,bTR)*
     &           (3*CB*dup18*MW*Mf(bTR,2)*UASfC(All5,5,bTR)*
     &              ZNeu(Neu5,3) + 
     &             CB2*MW2*UASfC(All5,2,bTR)*
     &              (-3*dup18*ZNeu(Neu5,2) + 
     &                ZNeu(Neu5,1)*
     &                 (SW2*ZNeuC(Neu5,1) - 3*CW*SW*ZNeuC(Neu5,2)))
     &             ) + 3*Mf(bTR,3)*UASf(All5,6,bTR)*
     &           (CB*MW*UASfC(All5,2,bTR)*
     &              (CW*SW*ZNeu(Neu5,1) - 3*CW2*ZNeu(Neu5,2)) + 
     &             3*CW2*Mf(bTR,2)*UASfC(All5,5,bTR)*ZNeu(Neu5,3))*
     &           ZNeuC(Neu5,3)))/
     &      (CB2*CW2**2*MZ2*CKM(3,3)*CKMC(3,2)*
     &        (MASf2(All5,bTR) - MNeu2(Neu5))**2)

	ENDLOOP(All5)
	ENDLOOP(Neu5)

	LOOP(Neu6, 1,4,1)
	LOOP(Neu5, 1,4,1)
	LOOP(All5, 1,6,1)

	dup19 = CW*SW*ZNeu(Neu6,1) - 3*CW2*ZNeu(Neu6,2)

	dup20 = CW*SW*ZNeuC(Neu5,1) - 3*CW2*ZNeuC(Neu5,2)

	dup21 = SW2*ZNeuC(Neu5,1) - 3*CW*SW*ZNeuC(Neu5,2)

        dup22 = ZNeu(Neu5,3)*ZNeuC(Neu6,3) - 
     &    ZNeu(Neu5,4)*ZNeuC(Neu6,4)

        dup23 = CB2*MW2*UASfC(All5,2,bTR)*
     &     (dup21*ZNeu(Neu6,1) - 3*dup20*ZNeu(Neu6,2)) + 
     &    3*CB*dup20*MW*Mf(bTR,2)*UASfC(All5,5,bTR)*ZNeu(Neu6,3)

        CALNeu = CALNeu + 
     &    1/144.D0*(dup22*(-1 + 
     &           4*C00z(MASf2(All5,bTR),MNeu2(Neu5),MNeu2(Neu6)))*
     &         (dup23*UASf(All5,3,bTR) + 
     &           3*Mf(bTR,3)*UASf(All5,6,bTR)*
     &            (CB*dup19*MW*UASfC(All5,2,bTR) + 
     &              3*CW2*Mf(bTR,2)*UASfC(All5,5,bTR)*ZNeu(Neu6,3))
     &             *ZNeuC(Neu5,3)) + 
     &        2*C0z(MASf2(All5,bTR),MNeu2(Neu5),MNeu2(Neu6))*
     &         MNeu(Neu5)*(-(dup22*MB*
     &              (3*Mf(bTR,3)*UASf(All5,3,bTR)*ZNeu(Neu5,3)*
     &                 (CB*dup19*MW*UASfC(All5,2,bTR) + 
     &                   3*CW2*Mf(bTR,2)*UASfC(All5,5,bTR)*
     &                    ZNeu(Neu6,3)) + 
     &                2*UASf(All5,6,bTR)*ZNeu(Neu5,1)*
     &                 (CB2*MW2*UASfC(All5,2,bTR)*
     &                    (SW2*ZNeu(Neu6,1) - 3*CW*SW*ZNeu(Neu6,2))
     &                     + 3*CB*CW*MW*SW*Mf(bTR,2)*
     &                    UASfC(All5,5,bTR)*ZNeu(Neu6,3)))) + 
     &           MNeu(Neu6)*
     &            (dup23*UASf(All5,3,bTR) + 
     &              3*Mf(bTR,3)*UASf(All5,6,bTR)*
     &               (CB*dup19*MW*UASfC(All5,2,bTR) + 
     &                 3*CW2*Mf(bTR,2)*UASfC(All5,5,bTR)*
     &                  ZNeu(Neu6,3))*ZNeuC(Neu5,3))*
     &            (ZNeu(Neu6,3)*ZNeuC(Neu5,3) - 
     &              ZNeu(Neu6,4)*ZNeuC(Neu5,4))))/
     &      (CB2*CW2**2*MZ2*CKM(3,3)*CKMC(3,2))

	ENDLOOP(All5)
	ENDLOOP(Neu5)
	ENDLOOP(Neu6)

	LOOP(Neu5, 1,4,1)
	LOOP(All6, 1,6,1)
	LOOP(All5, 1,6,1)

	tmp2 = C00z(MASf2(All5,bTR),MASf2(All6,bTR),MNeu2(Neu5))

	LOOP(Ind1, 1,3,1)

	dup24 = CW*SW*ZNeuC(Neu5,1) - 3*CW2*ZNeuC(Neu5,2)

        CALNeu = CALNeu + 
     &    1/108.D0*(tmp2*(UASf(All6,Ind1,bTR)*
     &           (UASfC(All5,Ind1,bTR) + 
     &             2*CW2*UASfC(All5,Ind1,bTR)) - 
     &          2*SW2*UASf(All6,3 + Ind1,bTR)*
     &           UASfC(All5,3 + Ind1,bTR))*
     &        (UASf(All5,3,bTR)*
     &           (3*CB*dup24*MW*Mf(bTR,2)*UASfC(All6,5,bTR)*
     &              ZNeu(Neu5,3) + 
     &             CB2*MW2*UASfC(All6,2,bTR)*
     &              (-3*dup24*ZNeu(Neu5,2) + 
     &                ZNeu(Neu5,1)*
     &                 (SW2*ZNeuC(Neu5,1) - 3*CW*SW*ZNeuC(Neu5,2)))
     &             ) + 3*Mf(bTR,3)*UASf(All5,6,bTR)*
     &           (CB*MW*UASfC(All6,2,bTR)*
     &              (CW*SW*ZNeu(Neu5,1) - 3*CW2*ZNeu(Neu5,2)) + 
     &             3*CW2*Mf(bTR,2)*UASfC(All6,5,bTR)*ZNeu(Neu5,3))*
     &           ZNeuC(Neu5,3)))/
     &      (CB2*CW2**2*MZ2*CKM(3,3)*CKMC(3,2))

	ENDLOOP(Ind1)

	ENDLOOP(All5)
	ENDLOOP(All6)
	ENDLOOP(Neu5)

	LOOP(Sfe5, 1,2,1)
	LOOP(Neu6, 1,4,1)
	LOOP(Neu5, 1,4,1)
	LOOP(All5, 1,6,1)

	dup25 = CW*SW*ZNeu(Neu6,1) - 3*CW2*ZNeu(Neu6,2)

	dup26 = CW*SW*ZNeu(Neu6,1) + CW2*ZNeu(Neu6,2)

	dup27 = CW*SW*ZNeuC(Neu5,1) - 3*CW2*ZNeuC(Neu5,2)

	dup28 = CW*SW*ZNeuC(Neu5,1) + CW2*ZNeuC(Neu5,2)

	dup29 = SW2*ZNeuC(Neu5,1) - 3*CW*SW*ZNeuC(Neu5,2)

	dup30 = SW2*ZNeuC(Neu5,1) + CW*SW*ZNeuC(Neu5,2)

	dup31 = CW*SW*ZNeuC(Neu6,1) + CW2*ZNeuC(Neu6,2)

	dup32 = SW2*ZNeuC(Neu6,1) + CW*SW*ZNeuC(Neu6,2)

        dup33 = CW2*ZNeu(Neu5,2)*
     &     (dup32*ZNeu(Neu6,1) - 3*dup31*ZNeu(Neu6,2)) + 
     &    ZNeu(Neu5,1)*(dup31*SW2*ZNeu(Neu6,1) - 
     &       3*CW2*dup32*ZNeu(Neu6,2))

        dup34 = USf(Sfe5,1,2,2)*USfC(Sfe5,1,2,2)*
     &     (dup32*ZNeu(Neu5,1) + dup31*ZNeu(Neu5,2)) + 
     &    4*SW2*USf(Sfe5,2,2,2)*USfC(Sfe5,2,2,2)*ZNeu(Neu5,1)*
     &     ZNeuC(Neu6,1)

        tmp3 = -2*D00z(MASf2(All5,bTR),MNeu2(Neu5),MNeu2(Neu6),
     &     MSf2(Sfe5,2,2))*
     &    (3*Mf(bTR,3)*UASf(All5,6,bTR)*ZNeuC(Neu5,3)*
     &       (3*CW2*dup34*Mf(bTR,2)*UASfC(All5,5,bTR)*
     &          ZNeu(Neu6,3) + 
     &         CB*MW*UASfC(All5,2,bTR)*
     &          (dup33*USf(Sfe5,1,2,2)*USfC(Sfe5,1,2,2) + 
     &            4*dup25*SW2*USf(Sfe5,2,2,2)*USfC(Sfe5,2,2,2)*
     &             ZNeu(Neu5,1)*ZNeuC(Neu6,1))) + 
     &      UASf(All5,3,bTR)*
     &       (3*CB*MW*Mf(bTR,2)*UASfC(All5,5,bTR)*ZNeu(Neu6,3)*
     &          (USf(Sfe5,1,2,2)*USfC(Sfe5,1,2,2)*
     &             (CW2*ZNeu(Neu5,2)*
     &                (dup32*ZNeuC(Neu5,1) - 3*dup31*ZNeuC(Neu5,2))
     &                 + ZNeu(Neu5,1)*
     &                (dup31*SW2*ZNeuC(Neu5,1) - 
     &                  3*CW2*dup32*ZNeuC(Neu5,2))) + 
     &            4*dup27*SW2*USf(Sfe5,2,2,2)*USfC(Sfe5,2,2,2)*
     &             ZNeu(Neu5,1)*ZNeuC(Neu6,1)) + 
     &         CB2*MW2*UASfC(All5,2,bTR)*
     &          (USf(Sfe5,1,2,2)*USfC(Sfe5,1,2,2)*
     &             (ZNeu(Neu5,2)*
     &                (3*CW2*ZNeu(Neu6,2)*
     &                   (-(dup32*ZNeuC(Neu5,1)) + 
     &                     3*dup31*ZNeuC(Neu5,2)) + 
     &                  ZNeu(Neu6,1)*
     &                   (dup31*SW2*ZNeuC(Neu5,1) - 
     &                     3*CW2*dup32*ZNeuC(Neu5,2))) + 
     &               ZNeu(Neu5,1)*
     &                (SW2*ZNeu(Neu6,1)*
     &                   (dup32*ZNeuC(Neu5,1) - 
     &                     3*dup31*ZNeuC(Neu5,2)) - 
     &                  3*ZNeu(Neu6,2)*
     &                   (dup31*SW2*ZNeuC(Neu5,1) - 
     &                     3*CW2*dup32*ZNeuC(Neu5,2)))) + 
     &            4*SW2*USf(Sfe5,2,2,2)*USfC(Sfe5,2,2,2)*
     &             ZNeu(Neu5,1)*
     &             (dup29*ZNeu(Neu6,1) - 3*dup27*ZNeu(Neu6,2))*
     &             ZNeuC(Neu6,1))))

        tmp3 = tmp3 + D0z(MASf2(All5,bTR),MNeu2(Neu5),
     &      MNeu2(Neu6),MSf2(Sfe5,2,2))*MNeu(Neu5)*
     &     (-(MNeu(Neu6)*(UASf(All5,3,bTR)*
     &             (3*CB*MW*Mf(bTR,2)*UASfC(All5,5,bTR)*
     &                ZNeu(Neu6,3)*
     &                (4*dup27*SW2*USf(Sfe5,2,2,2)*
     &                   USfC(Sfe5,2,2,2)*ZNeu(Neu6,1)*
     &                   ZNeuC(Neu5,1) + 
     &                  USf(Sfe5,1,2,2)*USfC(Sfe5,1,2,2)*
     &                   (-2*CW2*
     &                      (SW2*ZNeu(Neu6,1) + 
     &                      CW*SW*ZNeu(Neu6,2))*ZNeuC(Neu5,1)*
     &                      ZNeuC(Neu5,2) + 
     &                     dup26*
     &                      (SW2*ZNeuC(Neu5,1)**2 - 
     &                       3*CW2*ZNeuC(Neu5,2)**2))) + 
     &               CB2*MW2*UASfC(All5,2,bTR)*
     &                (4*SW2*USf(Sfe5,2,2,2)*USfC(Sfe5,2,2,2)*
     &                   ZNeu(Neu6,1)*
     &                   (dup29*ZNeu(Neu6,1) - 
     &                     3*dup27*ZNeu(Neu6,2))*ZNeuC(Neu5,1) + 
     &                  USf(Sfe5,1,2,2)*USfC(Sfe5,1,2,2)*
     &                   (SW2*ZNeu(Neu6,1)**2*
     &                      (SW2*ZNeuC(Neu5,1)**2 - 
     &                       2*CW*SW*ZNeuC(Neu5,1)*ZNeuC(Neu5,2) - 
     &                       3*CW2*ZNeuC(Neu5,2)**2) + 
     &                     3*CW2*ZNeu(Neu6,2)**2*
     &                      (-(SW2*ZNeuC(Neu5,1)**2) + 
     &                       2*CW*SW*ZNeuC(Neu5,1)*ZNeuC(Neu5,2) + 
     &                       3*CW2*ZNeuC(Neu5,2)**2) + 
     &                     2*ZNeu(Neu6,1)*ZNeu(Neu6,2)*
     &                      (2*CW2*SW2*ZNeuC(Neu5,1)*
     &                      ZNeuC(Neu5,2) + 
     &                       CW*SW*
     &                       (-(SW2*ZNeuC(Neu5,1)**2) + 
     &                       3*CW2*ZNeuC(Neu5,2)**2))))) + 
     &            3*Mf(bTR,3)*UASf(All5,6,bTR)*
     &             (3*CW2*Mf(bTR,2)*UASfC(All5,5,bTR)*ZNeu(Neu6,3)*
     &                (USf(Sfe5,1,2,2)*USfC(Sfe5,1,2,2)*
     &                   (dup30*ZNeu(Neu6,1) + dup28*ZNeu(Neu6,2))+
     &                    4*SW2*USf(Sfe5,2,2,2)*USfC(Sfe5,2,2,2)*
     &                   ZNeu(Neu6,1)*ZNeuC(Neu5,1)) + 
     &               CB*MW*UASfC(All5,2,bTR)*
     &                (USf(Sfe5,1,2,2)*USfC(Sfe5,1,2,2)*
     &                   (dup28*SW2*ZNeu(Neu6,1)**2 + 
     &                     CW2*
     &                      (-2*dup30*ZNeu(Neu6,1)*ZNeu(Neu6,2) - 
     &                       3*dup28*ZNeu(Neu6,2)**2)) + 
     &                  4*dup25*SW2*USf(Sfe5,2,2,2)*
     &                   USfC(Sfe5,2,2,2)*ZNeu(Neu6,1)*
     &                   ZNeuC(Neu5,1)))*ZNeuC(Neu5,3))) + 
     &       MB*(3*Mf(bTR,3)*UASf(All5,3,bTR)*ZNeu(Neu5,3)*
     &           (3*CW2*dup34*Mf(bTR,2)*UASfC(All5,5,bTR)*
     &              ZNeu(Neu6,3) + 
     &             CB*MW*UASfC(All5,2,bTR)*
     &              (dup33*USf(Sfe5,1,2,2)*USfC(Sfe5,1,2,2) + 
     &                4*dup25*SW2*USf(Sfe5,2,2,2)*USfC(Sfe5,2,2,2)*
     &                 ZNeu(Neu5,1)*ZNeuC(Neu6,1))) + 
     &          2*UASf(All5,6,bTR)*ZNeu(Neu5,1)*
     &           (3*CB*MW*Mf(bTR,2)*UASfC(All5,5,bTR)*ZNeu(Neu6,3)*
     &              (USf(Sfe5,1,2,2)*USfC(Sfe5,1,2,2)*
     &                 (dup31*SW2*ZNeu(Neu5,1) + 
     &                   CW2*dup32*ZNeu(Neu5,2)) + 
     &                4*CW*SW*SW2*USf(Sfe5,2,2,2)*USfC(Sfe5,2,2,2)*
     &                 ZNeu(Neu5,1)*ZNeuC(Neu6,1)) + 
     &             CB2*MW2*UASfC(All5,2,bTR)*
     &              (USf(Sfe5,1,2,2)*USfC(Sfe5,1,2,2)*
     &                 (SW2*ZNeu(Neu5,1)*
     &                    (dup32*ZNeu(Neu6,1) - 
     &                      3*dup31*ZNeu(Neu6,2)) + 
     &                   ZNeu(Neu5,2)*
     &                    (dup31*SW2*ZNeu(Neu6,1) - 
     &                      3*CW2*dup32*ZNeu(Neu6,2))) + 
     &                4*SW2*USf(Sfe5,2,2,2)*USfC(Sfe5,2,2,2)*
     &                 ZNeu(Neu5,1)*
     &                 (SW2*ZNeu(Neu6,1) - 3*CW*SW*ZNeu(Neu6,2))*
     &                 ZNeuC(Neu6,1)))))

        CALNeu = CALNeu + 
     &    1/72.D0*tmp3/(CB2*CW2**2*CKM(3,3)*CKMC(3,2))

	ENDLOOP(All5)
	ENDLOOP(Neu5)
	ENDLOOP(Neu6)
	ENDLOOP(Sfe5)

#ifdef DETAILED_DEBUG
	DCONST "CALNeu =", CALNeu ENDL
#endif

	CARNeu = 0

	LOOP(Neu5, 1,4,1)
	LOOP(All5, 1,6,1)

        CARNeu = CARNeu - 
     &    1/432.D0*(SW2*(MASf2(All5,bTR)**2 - 
     &          2*A0(MASf2(All5,bTR))*
     &           (MASf2(All5,bTR) - 2*MNeu2(Neu5)) - 
     &          MNeu2(Neu5)*(2*A0(MNeu2(Neu5)) + MNeu2(Neu5)))*
     &        (3*Mf(bTR,3)*UASf(All5,3,bTR)*ZNeu(Neu5,3)*
     &           (2*CB*CW*MW*SW*UASfC(All5,5,bTR)*ZNeuC(Neu5,1) + 
     &             3*CW2*Mf(bTR,2)*UASfC(All5,2,bTR)*ZNeuC(Neu5,3))
     &            + 2*UASf(All5,6,bTR)*ZNeu(Neu5,1)*
     &           (2*CB2*MW2*SW2*UASfC(All5,5,bTR)*ZNeuC(Neu5,1) + 
     &             3*CB*CW*MW*SW*Mf(bTR,2)*UASfC(All5,2,bTR)*
     &              ZNeuC(Neu5,3))))/
     &      (CB2*CW2**2*MZ2*CKM(3,3)*CKMC(3,2)*
     &        (MASf2(All5,bTR) - MNeu2(Neu5))**2)

	ENDLOOP(All5)
	ENDLOOP(Neu5)

	LOOP(Neu6, 1,4,1)
	LOOP(Neu5, 1,4,1)
	LOOP(All5, 1,6,1)

        dup35 = -(ZNeu(Neu6,3)*ZNeuC(Neu5,3)) + 
     &    ZNeu(Neu6,4)*ZNeuC(Neu5,4)

        dup36 = ZNeu(Neu5,3)*ZNeuC(Neu6,3) - 
     &    ZNeu(Neu5,4)*ZNeuC(Neu6,4)

	dup37 = 1 - 4*C00z(MASf2(All5,bTR),MNeu2(Neu5),MNeu2(Neu6))

	dup38 = -1 + 4*C00z(MASf2(All5,bTR),MNeu2(Neu5),MNeu2(Neu6))

        CARNeu = CARNeu - 
     &    1/144.D0*(2*UASf(All5,6,bTR)*
     &         (2*UASfC(All5,5,bTR)*
     &            (3*CB*CW*dup35*MB*MW*SW*
     &               C0z(MASf2(All5,bTR),MNeu2(Neu5),MNeu2(Neu6))*
     &               Mf(bTR,3)*MNeu(Neu5)*ZNeuC(Neu5,3) + 
     &              CB2*MW2*SW2*ZNeu(Neu5,1)*
     &               (2*dup36*
     &                  C0z(MASf2(All5,bTR),MNeu2(Neu5),
     &                   MNeu2(Neu6))*MNeu(Neu5)*MNeu(Neu6) + 
     &                 dup38*ZNeu(Neu6,3)*ZNeuC(Neu5,3) + 
     &                 dup37*ZNeu(Neu6,4)*ZNeuC(Neu5,4)))*
     &            ZNeuC(Neu6,1) + 
     &           3*Mf(bTR,2)*UASfC(All5,2,bTR)*
     &            (3*CW2*dup35*MB*
     &               C0z(MASf2(All5,bTR),MNeu2(Neu5),MNeu2(Neu6))*
     &               Mf(bTR,3)*MNeu(Neu5)*ZNeuC(Neu5,3) + 
     &              CB*CW*MW*SW*ZNeu(Neu5,1)*
     &               (2*dup36*
     &                  C0z(MASf2(All5,bTR),MNeu2(Neu5),
     &                   MNeu2(Neu6))*MNeu(Neu5)*MNeu(Neu6) + 
     &                 dup38*ZNeu(Neu6,3)*ZNeuC(Neu5,3) + 
     &                 dup37*ZNeu(Neu6,4)*ZNeuC(Neu5,4)))*
     &            ZNeuC(Neu6,3)) + 
     &        UASf(All5,3,bTR)*
     &         (3*Mf(bTR,3)*ZNeu(Neu5,3)*
     &            (2*dup36*
     &               C0z(MASf2(All5,bTR),MNeu2(Neu5),MNeu2(Neu6))*
     &               MNeu(Neu5)*MNeu(Neu6) + 
     &              dup38*ZNeu(Neu6,3)*ZNeuC(Neu5,3) + 
     &              dup37*ZNeu(Neu6,4)*ZNeuC(Neu5,4))*
     &            (2*CB*CW*MW*SW*UASfC(All5,5,bTR)*ZNeuC(Neu6,1) + 
     &              3*CW2*Mf(bTR,2)*UASfC(All5,2,bTR)*ZNeuC(Neu6,3)
     &              ) - 2*MB*
     &            C0z(MASf2(All5,bTR),MNeu2(Neu5),MNeu2(Neu6))*
     &            MNeu(Neu5)*
     &            (ZNeu(Neu6,3)*ZNeuC(Neu5,3) - 
     &              ZNeu(Neu6,4)*ZNeuC(Neu5,4))*
     &            (2*CB2*MW2*UASfC(All5,5,bTR)*
     &               (SW2*ZNeuC(Neu5,1) - 3*CW*SW*ZNeuC(Neu5,2))*
     &               ZNeuC(Neu6,1) + 
     &              3*CB*MW*Mf(bTR,2)*UASfC(All5,2,bTR)*
     &               (CW*SW*ZNeuC(Neu5,1) - 3*CW2*ZNeuC(Neu5,2))*
     &               ZNeuC(Neu6,3))))/
     &      (CB2*CW2**2*MZ2*CKM(3,3)*CKMC(3,2))

	ENDLOOP(All5)
	ENDLOOP(Neu5)
	ENDLOOP(Neu6)

	LOOP(Neu5, 1,4,1)
	LOOP(All6, 1,6,1)
	LOOP(All5, 1,6,1)

	tmp4 = C00z(MASf2(All5,bTR),MASf2(All6,bTR),MNeu2(Neu5))

	LOOP(Ind1, 1,3,1)

        CARNeu = CARNeu + 
     &    1/108.D0*(tmp4*(UASf(All6,Ind1,bTR)*
     &           (UASfC(All5,Ind1,bTR) + 
     &             2*CW2*UASfC(All5,Ind1,bTR)) - 
     &          2*SW2*UASf(All6,3 + Ind1,bTR)*
     &           UASfC(All5,3 + Ind1,bTR))*
     &        (3*Mf(bTR,3)*UASf(All5,3,bTR)*ZNeu(Neu5,3)*
     &           (2*CB*CW*MW*SW*UASfC(All6,5,bTR)*ZNeuC(Neu5,1) + 
     &             3*CW2*Mf(bTR,2)*UASfC(All6,2,bTR)*ZNeuC(Neu5,3))
     &            + 2*UASf(All5,6,bTR)*ZNeu(Neu5,1)*
     &           (2*CB2*MW2*SW2*UASfC(All6,5,bTR)*ZNeuC(Neu5,1) + 
     &             3*CB*CW*MW*SW*Mf(bTR,2)*UASfC(All6,2,bTR)*
     &              ZNeuC(Neu5,3))))/
     &      (CB2*CW2**2*MZ2*CKM(3,3)*CKMC(3,2))

	ENDLOOP(Ind1)

	ENDLOOP(All5)
	ENDLOOP(All6)
	ENDLOOP(Neu5)

	LOOP(Sfe5, 1,2,1)
	LOOP(Neu6, 1,4,1)
	LOOP(Neu5, 1,4,1)
	LOOP(All5, 1,6,1)

	dup39 = CW*SW*ZNeuC(Neu5,1) + CW2*ZNeuC(Neu5,2)

	dup40 = SW2*ZNeuC(Neu5,1) + CW*SW*ZNeuC(Neu5,2)

	dup41 = CW*SW*ZNeuC(Neu6,1) + CW2*ZNeuC(Neu6,2)

	dup42 = SW2*ZNeuC(Neu6,1) + CW*SW*ZNeuC(Neu6,2)

        dup43 = -(SW2*ZNeuC(Neu5,1)**2) + 
     &    2*CW*SW*ZNeuC(Neu5,1)*ZNeuC(Neu5,2) + 
     &    3*CW2*ZNeuC(Neu5,2)**2

        dup44 = 2*CW2*SW2*ZNeuC(Neu5,1)*ZNeuC(Neu5,2) + 
     &    CW*SW*(-(SW2*ZNeuC(Neu5,1)**2) + 3*CW2*ZNeuC(Neu5,2)**2)

        dup45 = USf(Sfe5,1,2,2)*USfC(Sfe5,1,2,2)*
     &     (dup40*ZNeu(Neu6,1) + dup39*ZNeu(Neu6,2)) + 
     &    4*SW2*USf(Sfe5,2,2,2)*USfC(Sfe5,2,2,2)*ZNeu(Neu6,1)*
     &     ZNeuC(Neu5,1)

        dup46 = USf(Sfe5,1,2,2)*USfC(Sfe5,1,2,2)*
     &     (dup39*SW2*ZNeu(Neu6,1) + CW2*dup40*ZNeu(Neu6,2)) + 
     &    4*CW*SW*SW2*USf(Sfe5,2,2,2)*USfC(Sfe5,2,2,2)*
     &     ZNeu(Neu6,1)*ZNeuC(Neu5,1)

        CARNeu = CARNeu + 
     &    1/72.D0*(2*D00z(MASf2(All5,bTR),MNeu2(Neu5),MNeu2(Neu6),
     &          MSf2(Sfe5,2,2))*
     &         (3*Mf(bTR,3)*UASf(All5,3,bTR)*ZNeu(Neu5,3)*
     &            (2*CB*dup46*MW*UASfC(All5,5,bTR)*ZNeuC(Neu6,1) + 
     &              3*CW2*dup45*Mf(bTR,2)*UASfC(All5,2,bTR)*
     &               ZNeuC(Neu6,3)) + 
     &           2*UASf(All5,6,bTR)*ZNeu(Neu5,1)*
     &            (2*CB2*dup45*MW2*SW2*UASfC(All5,5,bTR)*
     &               ZNeuC(Neu6,1) + 
     &              3*CB*dup46*MW*Mf(bTR,2)*UASfC(All5,2,bTR)*
     &               ZNeuC(Neu6,3))) + 
     &        D0z(MASf2(All5,bTR),MNeu2(Neu5),MNeu2(Neu6),
     &          MSf2(Sfe5,2,2))*MNeu(Neu5)*
     &         (UASf(All5,3,bTR)*
     &            (2*UASfC(All5,5,bTR)*ZNeuC(Neu6,1)*
     &               (USf(Sfe5,1,2,2)*USfC(Sfe5,1,2,2)*
     &                  (3*CB*MW*Mf(bTR,3)*MNeu(Neu6)*
     &                     (dup41*SW2*ZNeu(Neu5,1) + 
     &                       CW2*dup42*ZNeu(Neu5,2))*ZNeu(Neu5,3)+
     &                      CB2*MB*MW2*
     &                     (dup43*SW2*ZNeu(Neu6,1) + 
     &                       dup44*ZNeu(Neu6,2))) + 
     &                 4*SW2*USf(Sfe5,2,2,2)*USfC(Sfe5,2,2,2)*
     &                  (CB2*MB*MW2*ZNeu(Neu6,1)*ZNeuC(Neu5,1)*
     &                     (-(SW2*ZNeuC(Neu5,1)) + 
     &                       3*CW*SW*ZNeuC(Neu5,2)) + 
     &                    3*CB*CW*MW*SW*Mf(bTR,3)*MNeu(Neu6)*
     &                     ZNeu(Neu5,1)*ZNeu(Neu5,3)*ZNeuC(Neu6,1))
     &                 ) + 
     &              3*Mf(bTR,2)*UASfC(All5,2,bTR)*
     &               (USf(Sfe5,1,2,2)*USfC(Sfe5,1,2,2)*
     &                  (CB*dup44*MB*MW*ZNeu(Neu6,1) + 
     &                    CW2*
     &                     (3*Mf(bTR,3)*MNeu(Neu6)*
     &                       (dup42*ZNeu(Neu5,1) + 
     &                       dup41*ZNeu(Neu5,2))*ZNeu(Neu5,3) + 
     &                       CB*dup43*MB*MW*ZNeu(Neu6,2))) + 
     &                 4*SW2*USf(Sfe5,2,2,2)*USfC(Sfe5,2,2,2)*
     &                  (CB*MB*MW*ZNeu(Neu6,1)*ZNeuC(Neu5,1)*
     &                     (-(CW*SW*ZNeuC(Neu5,1)) + 
     &                       3*CW2*ZNeuC(Neu5,2)) + 
     &                    3*CW2*Mf(bTR,3)*MNeu(Neu6)*ZNeu(Neu5,1)*
     &                     ZNeu(Neu5,3)*ZNeuC(Neu6,1)))*
     &               ZNeuC(Neu6,3)) + 
     &           UASf(All5,6,bTR)*
     &            (-3*MB*Mf(bTR,3)*ZNeuC(Neu5,3)*
     &               (2*CB*dup46*MW*UASfC(All5,5,bTR)*
     &                  ZNeuC(Neu6,1) + 
     &                 3*CW2*dup45*Mf(bTR,2)*UASfC(All5,2,bTR)*
     &                  ZNeuC(Neu6,3)) + 
     &              2*MNeu(Neu6)*ZNeu(Neu5,1)*
     &               (2*CB2*MW2*SW2*UASfC(All5,5,bTR)*
     &                  ZNeuC(Neu6,1)*
     &                  (USf(Sfe5,1,2,2)*USfC(Sfe5,1,2,2)*
     &                     (dup42*ZNeu(Neu5,1) + 
     &                       dup41*ZNeu(Neu5,2)) + 
     &                    4*SW2*USf(Sfe5,2,2,2)*USfC(Sfe5,2,2,2)*
     &                     ZNeu(Neu5,1)*ZNeuC(Neu6,1)) + 
     &                 3*CB*MW*Mf(bTR,2)*UASfC(All5,2,bTR)*
     &                  (USf(Sfe5,1,2,2)*USfC(Sfe5,1,2,2)*
     &                     (dup41*SW2*ZNeu(Neu5,1) + 
     &                       CW2*dup42*ZNeu(Neu5,2)) + 
     &                    4*CW*SW*SW2*USf(Sfe5,2,2,2)*
     &                     USfC(Sfe5,2,2,2)*ZNeu(Neu5,1)*
     &                     ZNeuC(Neu6,1))*ZNeuC(Neu6,3)))))/
     &      (CB2*CW2**2*CKM(3,3)*CKMC(3,2))

	ENDLOOP(All5)
	ENDLOOP(Neu5)
	ENDLOOP(Neu6)
	ENDLOOP(Sfe5)

#ifdef DETAILED_DEBUG
	DCONST "CARNeu =", CARNeu ENDL
#endif

