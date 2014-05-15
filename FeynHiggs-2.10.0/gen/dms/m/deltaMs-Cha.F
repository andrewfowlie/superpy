
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

