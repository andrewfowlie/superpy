
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

