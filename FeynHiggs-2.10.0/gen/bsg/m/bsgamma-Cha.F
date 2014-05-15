
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

