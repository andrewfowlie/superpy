
	C7LSM = 0

	LOOP(Gen4, 1,3,1)

        dup1 = 2*(MB2 - MW2)*MW2 + 
     &    Mf2(3,Gen4)*(-MB2 + MW2 + Mf2(3,Gen4))

        C7LSM = C7LSM - 
     &    1/12.D0*(CKM(Gen4,3)*CKMC(Gen4,2)*
     &        (2*dup1*A0(Mf2(3,Gen4)) - 
     &          2*A0(MW2)*(2*(4*MB2 - MW2)*MW2 + 
     &             Mf2(3,Gen4)*(-7*MB2 + MW2 + Mf2(3,Gen4))) + 
     &          (MW2 - Mf2(3,Gen4))*
     &           (2*dup1*B0(MB2,MW2,Mf2(3,Gen4)) + 
     &             MB2*(MW2*
     &                 (2 - 
     &                   6*C0(MB2,0.D0,0.D0,MW2,Mf2(3,Gen4),MW2)*
     &                    (-2*MB2 + 2*MW2 + Mf2(3,Gen4))) - 
     &                Mf2(3,Gen4)*
     &                 (5 + 
     &                   4*
     &                    C0(MB2,0.D0,0.D0,MW2,Mf2(3,Gen4),
     &                     Mf2(3,Gen4))*
     &                    (-MB2 + 2*MW2 + Mf2(3,Gen4)))))))/
     &      (MB2**2*CKM(3,3)*CKMC(3,2)*(MW2 - Mf2(3,Gen4)))

	ENDLOOP(Gen4)

#ifdef DETAILED_DEBUG
	DCONST "C7LSM =", C7LSM ENDL
#endif


	C7LHp = 0

	LOOP(Gen4, 1,3,1)

	dup1 = MHp2 - MB*TB2*Mf(bTR,3) - Mf2(3,Gen4)

        C7LHp = C7LHp + 
     &    1/12.D0*(CKM(Gen4,3)*CKMC(Gen4,2)*Mf2(3,Gen4)*
     &        (2*(dup1*A0(MHp2) + 
     &             A0(Mf2(3,Gen4))*
     &              (-MHp2 + MB*TB2*Mf(bTR,3) + Mf2(3,Gen4))) - 
     &          (MHp2 - Mf2(3,Gen4))*
     &           (2*dup1*B0(MB2,MHp2,Mf2(3,Gen4)) + 
     &             MB2*(5 + 
     &                6*MHp2*C0(MB2,0.D0,0.D0,MHp2,Mf2(3,Gen4),MHp2) + 
     &                4*C0(MB2,0.D0,0.D0,MHp2,Mf2(3,Gen4),Mf2(3,Gen4))*
     &                 (MB*TB2*Mf(bTR,3) + Mf2(3,Gen4))))))/
     &      (MB2**2*TB2*CKM(3,3)*CKMC(3,2)*(-MHp2 + Mf2(3,Gen4)))

	ENDLOOP(Gen4)

#ifdef DETAILED_DEBUG
	DCONST "C7LHp =", C7LHp ENDL
#endif


	C7LCha = 0

	LOOP(All4, 1,6,1)

	tmp5 = A0(MASf2(All4,3))

	LOOP(Cha4, 1,2,1)

	tmp1 = A0(MCha2(Cha4))

	tmp2 = B0(MB2,MASf2(All4,3),MCha2(Cha4))

	tmp3 = C0(MB2,0.D0,0.D0,MASf2(All4,3),MCha2(Cha4),MASf2(All4,3))

	tmp4 = C0(MB2,0.D0,0.D0,MASf2(All4,3),MCha2(Cha4),MCha2(Cha4))

	LOOP(Ind2, 1,3,1)
	LOOP(Ind1, 1,3,1)

        dup1 = -(sqrt2*
     &       (Mf(3,Ind2)*UASfC(All4,3 + Ind2,3)*VCha(Cha4,2))) + 
     &    2*MW*SB*UASfC(All4,Ind2,3)*VCha(Cha4,1)

        dup2 = -(sqrt2*
     &       (Mf(3,Ind1)*UASf(All4,3 + Ind1,3)*VChaC(Cha4,2))) + 
     &    2*MW*SB*UASf(All4,Ind1,3)*VChaC(Cha4,1)

        C7LCha = C7LCha + 
     &    1/12.D0*((1/2.D0*(dup1*dup2*CKM(Ind1,3)*CKMC(Ind2,2)*
     &               (2*tmp1 - 2*tmp5 - 
     &                 MB2*
     &                  (5 + 4*tmp3*MASf2(All4,3) + 
     &                    6*tmp4*MASf2(All4,3)) + 
     &                 tmp2*(MB2 + MASf2(All4,3) - MCha2(Cha4)) - 
     &                 tmp2*(MB2 - MASf2(All4,3) + MCha2(Cha4))))/
     &             SB2 + sqrt2*
     &            (MB*CKM(Ind2,3)*CKMC(Ind1,2)*MCha(Cha4)*
     &               (-tmp2 + 
     &                 (-tmp1 + tmp5)/(MASf2(All4,3) - MCha2(Cha4))
     &                 )*Mf(bTR,3)*UASf(All4,Ind2,3)*UCha(Cha4,2)*
     &               (sqrt2*
     &                  (Mf(3,Ind1)*UASfC(All4,3 + Ind1,3)*
     &                    VCha(Cha4,2)) - 
     &                 2*MW*SB*UASfC(All4,Ind1,3)*VCha(Cha4,1)))/
     &             (CB*SB))/MB2 - 
     &        3*tmp4*((2*MB*CKM(Ind2,3)*CKMC(Ind1,2)*MCha(Cha4)*
     &              Mf(bTR,3)*UASf(All4,Ind2,3)*UCha(Cha4,2)*
     &              (-(sqrt2*
     &                   (MW*SB*UASfC(All4,Ind1,3)*VCha(Cha4,1)))+
     &                  Mf(3,Ind1)*UASfC(All4,3 + Ind1,3)*
     &                 VCha(Cha4,2)))/(CB*SB) + 
     &           (CKM(Ind1,3)*CKMC(Ind2,2)*
     &              (-MASf2(All4,3) + MCha2(Cha4))*
     &              (sqrt2*
     &                 (Mf(3,Ind2)*UASfC(All4,3 + Ind2,3)*
     &                   VCha(Cha4,2)) - 
     &                2*MW*SB*UASfC(All4,Ind2,3)*VCha(Cha4,1))*
     &              (sqrt2*
     &                 (Mf(3,Ind1)*UASf(All4,3 + Ind1,3)*
     &                   VChaC(Cha4,2)) - 
     &                2*MW*SB*UASf(All4,Ind1,3)*VChaC(Cha4,1)))/SB2
     &           ))/(MB2*CKM(3,3)*CKMC(3,2))

	ENDLOOP(Ind1)
	ENDLOOP(Ind2)

	ENDLOOP(Cha4)

	ENDLOOP(All4)

#ifdef DETAILED_DEBUG
	DCONST "C7LCha =", C7LCha ENDL
#endif


	C7LNeu = 0

	LOOP(All4, 1,6,1)

	tmp1 = A0(MASf2(All4,bTR))

	LOOP(Neu4, 1,4,1)

        dup1 = 2*CB*MW*SW*UASf(All4,6,bTR)*ZNeu(Neu4,1) + 
     &    3*CW*Mf(bTR,3)*UASf(All4,3,bTR)*ZNeu(Neu4,3)

        C7LNeu = C7LNeu + 
     &    1/54.D0*(MW2*(1/2.D0*
     &           ((MB2 - 2*tmp1 + 
     &               2*(A0(MNeu2(Neu4)) + 
     &                  MB2*
     &                   C0(MB2,0.D0,0.D0,MASf2(All4,bTR),MNeu2(Neu4),
     &                    MASf2(All4,bTR))*MASf2(All4,bTR) + 
     &                  B0(MB2,MASf2(All4,bTR),MNeu2(Neu4))*
     &                   (MASf2(All4,bTR) - MNeu2(Neu4))))*
     &             (CB*MW*UASf(All4,3,bTR)*
     &                (SW*ZNeuC(Neu4,1) - 3*CW*ZNeuC(Neu4,2)) + 
     &               3*CW*Mf(bTR,3)*UASf(All4,6,bTR)*ZNeuC(Neu4,3))
     &             ) + dup1*MB*MNeu(Neu4)*
     &           (-B0(MB2,MASf2(All4,bTR),MNeu2(Neu4)) + 
     &             (tmp1 - A0(MNeu2(Neu4)))/
     &              (MASf2(All4,bTR) - MNeu2(Neu4))))*
     &        UASfC(All4,2,bTR)*
     &        (SW*ZNeu(Neu4,1) - 3*CW*ZNeu(Neu4,2)))/
     &      (CB*CW2*MB2**2*MW*CKM(3,3)*CKMC(3,2))

	ENDLOOP(Neu4)

	ENDLOOP(All4)

#ifdef DETAILED_DEBUG
	DCONST "C7LNeu =", C7LNeu ENDL
#endif

	C7RNeu = 0

	LOOP(All4, 1,6,1)

	tmp2 = A0(MASf2(All4,bTR))

	LOOP(Neu4, 1,4,1)

        dup2 = 2*CB*MW*SW*UASf(All4,6,bTR)*ZNeu(Neu4,1) + 
     &    3*CW*Mf(bTR,3)*UASf(All4,3,bTR)*ZNeu(Neu4,3)

        C7RNeu = C7RNeu + 
     &    1/27.D0*(MW2*SW2*UASfC(All4,5,bTR)*ZNeuC(Neu4,1)*
     &        (1/2.D0*(dup2*(MB2 - 2*tmp2 + 
     &               2*(A0(MNeu2(Neu4)) + 
     &                  MB2*
     &                   C0(MB2,0.D0,0.D0,MASf2(All4,bTR),MNeu2(Neu4),
     &                    MASf2(All4,bTR))*MASf2(All4,bTR)) + 
     &               B0(MB2,MASf2(All4,bTR),MNeu2(Neu4))*
     &                (2*MASf2(All4,bTR) - 2*MNeu2(Neu4)))) + 
     &          MB*MNeu(Neu4)*
     &           (-B0(MB2,MASf2(All4,bTR),MNeu2(Neu4)) + 
     &             (tmp2 - A0(MNeu2(Neu4)))/
     &              (MASf2(All4,bTR) - MNeu2(Neu4)))*
     &           (CB*MW*UASf(All4,3,bTR)*
     &              (SW*ZNeuC(Neu4,1) - 3*CW*ZNeuC(Neu4,2)) + 
     &             3*CW*Mf(bTR,3)*UASf(All4,6,bTR)*ZNeuC(Neu4,3))))
     &       /(CB*CW2*MB2**2*MW*SW*CKM(3,3)*CKMC(3,2))

	ENDLOOP(Neu4)

	ENDLOOP(All4)

#ifdef DETAILED_DEBUG
	DCONST "C7RNeu =", C7RNeu ENDL
#endif


	C7LGlu = 0

	LOOP(All4, 1,6,1)

        C7LGlu = C7LGlu + 
     &    (2*Pi*sqrt2)/9.D0*(asMW*
     &        ((MGl2 - MASf2(All4,bTR))*
     &           (MB2 - 2*A0(MASf2(All4,bTR)) + 
     &             2*(A0(MGl2) + 
     &                MB2*C0(MB2,0.D0,0.D0,MGl2,MASf2(All4,bTR),
     &                  MASf2(All4,bTR))*MASf2(All4,bTR) + 
     &                B0(MB2,MGl2,MASf2(All4,bTR))*
     &                 (-MGl2 + MASf2(All4,bTR))))*UASf(All4,3,bTR)
     &            + 2*MB*M_3*
     &           (-A0(MGl2) + A0(MASf2(All4,bTR)) + 
     &             B0(MB2,MGl2,MASf2(All4,bTR))*
     &              (MGl2 - MASf2(All4,bTR)))*UASf(All4,6,bTR))*
     &        UASfC(All4,2,bTR))/
     &      (GF*MB2**2*CKM(3,3)*CKMC(3,2)*(MGl2 - MASf2(All4,bTR)))

	ENDLOOP(All4)

#ifdef DETAILED_DEBUG
	DCONST "C7LGlu =", C7LGlu ENDL
#endif

	C7RGlu = 0

	LOOP(All4, 1,6,1)

        C7RGlu = C7RGlu + 
     &    (2*Pi*sqrt2)/9.D0*(asMW*
     &        (2*MB*M_3C*
     &           (-A0(MGl2) + A0(MASf2(All4,bTR)) + 
     &             B0(MB2,MGl2,MASf2(All4,bTR))*
     &              (MGl2 - MASf2(All4,bTR)))*UASf(All4,3,bTR) + 
     &          (MGl2 - MASf2(All4,bTR))*
     &           (MB2 - 2*A0(MASf2(All4,bTR)) + 
     &             2*(A0(MGl2) + 
     &                MB2*C0(MB2,0.D0,0.D0,MGl2,MASf2(All4,bTR),
     &                  MASf2(All4,bTR))*MASf2(All4,bTR) + 
     &                B0(MB2,MGl2,MASf2(All4,bTR))*
     &                 (-MGl2 + MASf2(All4,bTR))))*UASf(All4,6,bTR)
     &          )*UASfC(All4,5,bTR))/
     &      (GF*MB2**2*CKM(3,3)*CKMC(3,2)*(MGl2 - MASf2(All4,bTR)))

	ENDLOOP(All4)

#ifdef DETAILED_DEBUG
	DCONST "C7RGlu =", C7RGlu ENDL
#endif

