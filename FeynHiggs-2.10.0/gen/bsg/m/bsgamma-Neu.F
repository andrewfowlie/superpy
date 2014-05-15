
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

