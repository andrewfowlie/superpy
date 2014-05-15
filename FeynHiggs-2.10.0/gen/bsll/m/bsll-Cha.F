
	CALCha = 0

	LOOP(Ind2, 1,3,1)
	LOOP(Ind1, 1,3,1)
	LOOP(Cha5, 1,2,1)
	LOOP(All5, 1,6,1)

        CALCha = CALCha - 
     &    1/96.D0*((1 + 2*CW2)*CKM(Ind1,3)*CKMC(Ind2,2)*
     &        (-MASf2(All5,3)**2 + 
     &          2*A0(MASf2(All5,3))*
     &           (MASf2(All5,3) - 2*MCha2(Cha5)) + 
     &          MCha2(Cha5)*(2*A0(MCha2(Cha5)) + MCha2(Cha5)))*
     &        (UASf(All5,Ind1,3)*
     &           (-(sqrt2*(MW*SB*Mf(3,Ind2)*UASfC(All5,3 + Ind2,3)*
     &                  VCha(Cha5,2))) + 
     &             2*MW2*SB2*UASfC(All5,Ind2,3)*VCha(Cha5,1))*
     &           VChaC(Cha5,1) + 
     &          Mf(3,Ind1)*UASf(All5,3 + Ind1,3)*
     &           (-(sqrt2*
     &                (MW*SB*UASfC(All5,Ind2,3)*VCha(Cha5,1))) + 
     &             Mf(3,Ind2)*UASfC(All5,3 + Ind2,3)*VCha(Cha5,2))*
     &           VChaC(Cha5,2)))/
     &      (CW2*MZ2*SB2*CKM(3,3)*CKMC(3,2)*
     &        (MASf2(All5,3) - MCha2(Cha5))**2)

	ENDLOOP(All5)
	ENDLOOP(Cha5)
	ENDLOOP(Ind1)
	ENDLOOP(Ind2)

	LOOP(Cha6, 1,2,1)
	LOOP(Cha5, 1,2,1)
	LOOP(All5, 1,6,1)

	tmp1 = C00z(MASf2(All5,3),MCha2(Cha5),MCha2(Cha6))

	tmp2 = C0z(MASf2(All5,3),MCha2(Cha5),MCha2(Cha6))

	tmp3 = D00z(MASf2(All5,3),MCha2(Cha5),MCha2(Cha6),MSf2(1,1,2))

	tmp4 = D0z(MASf2(All5,3),MCha2(Cha5),MCha2(Cha6),MSf2(1,1,2))

	LOOP(Ind2, 1,3,1)
	LOOP(Ind1, 1,3,1)

        dup1 = 2*SW2*Delta(Cha5,Cha6) - 
     &    2*VCha(Cha5,1)*VChaC(Cha6,1) - VCha(Cha5,2)*VChaC(Cha6,2)

        dup2 = sqrt2*(MW2*SB*UASfC(All5,Ind2,3)*
     &       VCha(Cha6,1)) - 
     &    MW*Mf(3,Ind2)*UASfC(All5,3 + Ind2,3)*VCha(Cha6,2)

        dup3 = sqrt2*(SB*Mf(3,Ind2)*UASfC(All5,3 + Ind2,3)*
     &       VCha(Cha6,2)) - 
     &    2*MW*SB2*UASfC(All5,Ind2,3)*VCha(Cha6,1)

        CALCha = CALCha + 
     &    1/(16.D0*sqrt2)*(MCha(Cha5)*
     &         ((-2*tmp2*(dup1*MB*SB2*CKM(Ind2,3)*CKMC(Ind1,2)*
     &                 Mf(bTR,3)*UASf(All5,Ind2,3)*UCha(Cha5,2)*
     &                 (-(sqrt2*
     &                      (Mf(3,Ind1)*UASfC(All5,3 + Ind1,3)*
     &                       VCha(Cha6,2))) + 
     &                   2*MW*SB*UASfC(All5,Ind1,3)*VCha(Cha6,1))-
     &                  CB*CKM(Ind1,3)*CKMC(Ind2,2)*MCha(Cha6)*
     &                 (2*SW2*Delta(Cha5,Cha6) - 
     &                   2*UCha(Cha6,1)*UChaC(Cha5,1) - 
     &                   UCha(Cha6,2)*UChaC(Cha5,2))*
     &                 (2*dup2*SB2*UASf(All5,Ind1,3)*
     &                    VChaC(Cha5,1) + 
     &                   dup3*Mf(3,Ind1)*UASf(All5,3 + Ind1,3)*
     &                    VChaC(Cha5,2))))/MZ2 + 
     &           4*CW2*MB*SB2*tmp4*CKM(Ind2,3)*CKMC(Ind1,2)*
     &            Mf(bTR,3)*UASf(All5,Ind2,3)*UCha(Cha5,2)*
     &            VCha(Cha5,1)*
     &            (sqrt2*(Mf(3,Ind1)*UASfC(All5,3 + Ind1,3)*
     &                 VCha(Cha6,2)) - 
     &              2*MW*SB*UASfC(All5,Ind1,3)*VCha(Cha6,1))*
     &            VChaC(Cha6,1)) + 
     &        (CB*CKM(Ind1,3)*CKMC(Ind2,2)*
     &           (2*dup2*SB2*UASf(All5,Ind1,3)*VChaC(Cha5,1) + 
     &             dup3*Mf(3,Ind1)*UASf(All5,3 + Ind1,3)*
     &              VChaC(Cha5,2))*
     &           (dup1*(-1 + SW2*(2 - 8*tmp1) + 4*tmp1) + 
     &             CW2*(dup1*(2 - 8*tmp1) - 
     &                8*MZ2*tmp3*VCha(Cha5,1)*VChaC(Cha6,1))))/MZ2)
     &       /(CB*CW2*SB*SB2*CKM(3,3)*CKMC(3,2))

	ENDLOOP(Ind1)
	ENDLOOP(Ind2)

	ENDLOOP(All5)
	ENDLOOP(Cha5)
	ENDLOOP(Cha6)

	LOOP(Cha5, 1,2,1)
	LOOP(All6, 1,6,1)
	LOOP(All5, 1,6,1)

	tmp5 = C00z(MASf2(All5,3),MASf2(All6,3),MCha2(Cha5))

	LOOP(Ind3, 1,3,1)
	LOOP(Ind2, 1,3,1)
	LOOP(Ind1, 1,3,1)

        CALCha = CALCha - 
     &    1/12.D0*(tmp5*CKM(Ind1,3)*CKMC(Ind2,2)*
     &        ((-1 + 4*CW2)*UASf(All6,Ind3,3)*UASfC(All5,Ind3,3) - 
     &          4*SW2*UASf(All6,3 + Ind3,3)*UASfC(All5,3 + Ind3,3))
     &         *(UASf(All5,Ind1,3)*
     &           (-(sqrt2*(MW*SB*Mf(3,Ind2)*UASfC(All6,3 + Ind2,3)*
     &                  VCha(Cha5,2))) + 
     &             2*MW2*SB2*UASfC(All6,Ind2,3)*VCha(Cha5,1))*
     &           VChaC(Cha5,1) + 
     &          Mf(3,Ind1)*UASf(All5,3 + Ind1,3)*
     &           (-(sqrt2*
     &                (MW*SB*UASfC(All6,Ind2,3)*VCha(Cha5,1))) + 
     &             Mf(3,Ind2)*UASfC(All6,3 + Ind2,3)*VCha(Cha5,2))*
     &           VChaC(Cha5,2)))/(CW2*MZ2*SB2*CKM(3,3)*CKMC(3,2))

	ENDLOOP(Ind1)
	ENDLOOP(Ind2)
	ENDLOOP(Ind3)

	ENDLOOP(All5)
	ENDLOOP(All6)
	ENDLOOP(Cha5)

#ifdef DETAILED_DEBUG
	DCONST "CALCha =", CALCha ENDL
#endif

	CARCha = 0

	LOOP(Ind2, 1,3,1)
	LOOP(Ind1, 1,3,1)
	LOOP(Cha5, 1,2,1)
	LOOP(All5, 1,6,1)

        CARCha = CARCha - 
     &    1/48.D0*(SW2*CKM(Ind1,3)*CKMC(Ind2,2)*
     &        (MASf2(All5,3)**2 - 
     &          2*A0(MASf2(All5,3))*
     &           (MASf2(All5,3) - 2*MCha2(Cha5)) - 
     &          MCha2(Cha5)*(2*A0(MCha2(Cha5)) + MCha2(Cha5)))*
     &        Mf(bTR,2)*Mf(bTR,3)*UASf(All5,Ind1,3)*
     &        UASfC(All5,Ind2,3)*UCha(Cha5,2)*UChaC(Cha5,2))/
     &      (CB2*CW2*MZ2*CKM(3,3)*CKMC(3,2)*
     &        (MASf2(All5,3) - MCha2(Cha5))**2)

	ENDLOOP(All5)
	ENDLOOP(Cha5)
	ENDLOOP(Ind1)
	ENDLOOP(Ind2)

	LOOP(Cha6, 1,2,1)
	LOOP(Cha5, 1,2,1)
	LOOP(All5, 1,6,1)

	tmp6 = C00z(MASf2(All5,3),MCha2(Cha5),MCha2(Cha6))

	tmp7 = C0z(MASf2(All5,3),MCha2(Cha5),MCha2(Cha6))

	tmp8 = D0z(MASf2(All5,3),MCha2(Cha5),MCha2(Cha6),MSf2(1,1,2))

	LOOP(Ind2, 1,3,1)
	LOOP(Ind1, 1,3,1)

	dup4 = 1 - 4*tmp6

        dup5 = -(sqrt2*
     &       (Mf(3,Ind1)*UASf(All5,3 + Ind1,3)*VChaC(Cha5,2))) + 
     &    2*MW*SB*UASf(All5,Ind1,3)*VChaC(Cha5,1)

        CARCha = CARCha + 
     &    1/(16.D0*sqrt2)*(CKM(Ind1,3)*CKMC(Ind2,2)*Mf(bTR,2)*
     &        UASfC(All5,Ind2,3)*
     &        (4*sqrt2*(CB*CW2*SB*tmp8*MCha(Cha5)*MCha(Cha6)*
     &             Mf(bTR,3)*UASf(All5,Ind1,3)*UCha(Cha5,2)*
     &             VCha(Cha5,1)*VChaC(Cha6,1)) - 
     &          (sqrt2*(CB*SB*Mf(bTR,3)*UASf(All5,Ind1,3)*
     &                UCha(Cha5,2)*
     &                (dup4*
     &                   (2*UCha(Cha6,1)*UChaC(Cha5,1) + 
     &                     UCha(Cha6,2)*UChaC(Cha5,2)) + 
     &                  2*tmp7*MCha(Cha5)*MCha(Cha6)*
     &                   (2*VCha(Cha5,1)*VChaC(Cha6,1) + 
     &                     VCha(Cha5,2)*VChaC(Cha6,2)))) + 
     &             2*SW2*Delta(Cha5,Cha6)*
     &              (sqrt2*
     &                 (CB*SB*
     &                   (-1 + 4*tmp6 - 
     &                     2*tmp7*MCha(Cha5)*MCha(Cha6))*Mf(bTR,3)*
     &                   UASf(All5,Ind1,3)*UCha(Cha5,2)) + 
     &                2*CB2*dup5*MB*tmp7*MCha(Cha5)) - 
     &             2*CB2*dup5*MB*tmp7*MCha(Cha5)*
     &              (2*UCha(Cha6,1)*UChaC(Cha5,1) + 
     &                UCha(Cha6,2)*UChaC(Cha5,2)))/MZ2)*
     &        UChaC(Cha6,2))/(CB*CB2*CW2*SB*CKM(3,3)*CKMC(3,2))

	ENDLOOP(Ind1)
	ENDLOOP(Ind2)

	ENDLOOP(All5)
	ENDLOOP(Cha5)
	ENDLOOP(Cha6)

	LOOP(Cha5, 1,2,1)
	LOOP(All6, 1,6,1)
	LOOP(All5, 1,6,1)

	tmp9 = C00z(MASf2(All5,3),MASf2(All6,3),MCha2(Cha5))

	LOOP(Ind3, 1,3,1)
	LOOP(Ind2, 1,3,1)
	LOOP(Ind1, 1,3,1)

        CARCha = CARCha - 
     &    1/12.D0*(tmp9*CKM(Ind2,3)*CKMC(Ind3,2)*Mf(bTR,2)*Mf(bTR,3)*
     &        UASf(All5,Ind2,3)*
     &        ((-1 + 4*CW2)*UASf(All6,Ind1,3)*UASfC(All5,Ind1,3) - 
     &          4*SW2*UASf(All6,3 + Ind1,3)*UASfC(All5,3 + Ind1,3))
     &         *UASfC(All6,Ind3,3)*UCha(Cha5,2)*UChaC(Cha5,2))/
     &      (CB2*CW2*MZ2*CKM(3,3)*CKMC(3,2))

	ENDLOOP(Ind1)
	ENDLOOP(Ind2)
	ENDLOOP(Ind3)

	ENDLOOP(All5)
	ENDLOOP(All6)
	ENDLOOP(Cha5)

#ifdef DETAILED_DEBUG
	DCONST "CARCha =", CARCha ENDL
#endif

