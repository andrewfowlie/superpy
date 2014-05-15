
	CALGlu = 0

	LOOP(All5, 1,6,1)

        CALGlu = CALGlu - 
     &    Pi/(18.D0*sqrt2)*(asMT*(1 + 2*CW2)*
     &        (MGl2*(MGl2 + 2*A0(MGl2) - 4*A0(MASf2(All5,bTR))) + 
     &          2*A0(MASf2(All5,bTR))*MASf2(All5,bTR) - 
     &          MASf2(All5,bTR)**2)*UASf(All5,3,bTR)*
     &        UASfC(All5,2,bTR))/
     &      (CW2*GF*MZ2*CKM(3,3)*CKMC(3,2)*
     &        (MGl2 - MASf2(All5,bTR))**2)

	ENDLOOP(All5)

	LOOP(All6, 1,6,1)
	LOOP(All5, 1,6,1)

	tmp1 = C00z(MGl2,MASf2(All5,bTR),MASf2(All6,bTR))

	LOOP(Ind1, 1,3,1)

        CALGlu = CALGlu + 
     &    (4*Pi)/(9.D0*sqrt2)*
     &     (asMT*tmp1*UASf(All5,3,bTR)*
     &        (UASf(All6,Ind1,bTR)*
     &           (UASfC(All5,Ind1,bTR) + 
     &             2*CW2*UASfC(All5,Ind1,bTR)) - 
     &          2*SW2*UASf(All6,3 + Ind1,bTR)*
     &           UASfC(All5,3 + Ind1,bTR))*UASfC(All6,2,bTR))/
     &      (CW2*GF*MZ2*CKM(3,3)*CKMC(3,2))

	ENDLOOP(Ind1)

	ENDLOOP(All5)
	ENDLOOP(All6)

#ifdef DETAILED_DEBUG
	DCONST "CALGlu =", CALGlu ENDL
#endif

	CARGlu = 0

	LOOP(All5, 1,6,1)

        CARGlu = CARGlu + 
     &    Pi/(9.D0*sqrt2)*(asMT*SW2*
     &        (MGl2*(MGl2 + 2*A0(MGl2) - 4*A0(MASf2(All5,bTR))) + 
     &          2*A0(MASf2(All5,bTR))*MASf2(All5,bTR) - 
     &          MASf2(All5,bTR)**2)*UASf(All5,6,bTR)*
     &        UASfC(All5,5,bTR))/
     &      (CW2*GF*MZ2*CKM(3,3)*CKMC(3,2)*
     &        (MGl2 - MASf2(All5,bTR))**2)

	ENDLOOP(All5)

	LOOP(All6, 1,6,1)
	LOOP(All5, 1,6,1)

	tmp2 = C00z(MGl2,MASf2(All5,bTR),MASf2(All6,bTR))

	LOOP(Ind1, 1,3,1)

        CARGlu = CARGlu + 
     &    (4*Pi)/(9.D0*sqrt2)*
     &     (asMT*tmp2*UASf(All5,6,bTR)*
     &        (UASf(All6,Ind1,bTR)*
     &           (UASfC(All5,Ind1,bTR) + 
     &             2*CW2*UASfC(All5,Ind1,bTR)) - 
     &          2*SW2*UASf(All6,3 + Ind1,bTR)*
     &           UASfC(All5,3 + Ind1,bTR))*UASfC(All6,5,bTR))/
     &      (CW2*GF*MZ2*CKM(3,3)*CKMC(3,2))

	ENDLOOP(Ind1)

	ENDLOOP(All5)
	ENDLOOP(All6)

#ifdef DETAILED_DEBUG
	DCONST "CARGlu =", CARGlu ENDL
#endif

