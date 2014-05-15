
	C8LGlu = 0

	LOOP(All4, 1,6,1)

        C8LGlu = C8LGlu + 
     &    Pi/(3.D0*sqrt2)*(asMW*
     &        ((MGl2 - MASf2(All4,bTR))*
     &           (-8*A0(MGl2) + 
     &             8*(A0(MASf2(All4,bTR)) + 
     &                B0(MB2,MGl2,MASf2(All4,bTR))*
     &                 (MGl2 - MASf2(All4,bTR))) + 
     &             MB2*(5 + 
     &                9*MGl2*
     &                 C0(MB2,0.D0,0.D0,MGl2,MASf2(All4,bTR),MGl2) + 
     &                C0(MB2,0.D0,0.D0,MGl2,MASf2(All4,bTR),
     &                  MASf2(All4,bTR))*MASf2(All4,bTR)))*
     &           UASf(All4,3,bTR) + 
     &          MB*M_3*(8*A0(MGl2) - 8*A0(MASf2(All4,bTR)) - 
     &             (8*B0(MB2,MGl2,MASf2(All4,bTR)) + 
     &                9*MB2*C0(MB2,0.D0,0.D0,MGl2,MASf2(All4,bTR),MGl2)
     &                )*(MGl2 - MASf2(All4,bTR)))*UASf(All4,6,bTR))
     &         *UASfC(All4,2,bTR))/
     &      (GF*MB2**2*CKM(3,3)*CKMC(3,2)*(MGl2 - MASf2(All4,bTR)))

	ENDLOOP(All4)

#ifdef DETAILED_DEBUG
	DCONST "C8LGlu =", C8LGlu ENDL
#endif

	C8RGlu = 0

	LOOP(All4, 1,6,1)

        C8RGlu = C8RGlu + 
     &    Pi/(3.D0*sqrt2)*(asMW*
     &        (MB*M_3C*(8*A0(MGl2) - 8*A0(MASf2(All4,bTR)) - 
     &             (8*B0(MB2,MGl2,MASf2(All4,bTR)) + 
     &                9*MB2*C0(MB2,0.D0,0.D0,MGl2,MASf2(All4,bTR),MGl2)
     &                )*(MGl2 - MASf2(All4,bTR)))*UASf(All4,3,bTR)+
     &            (MGl2 - MASf2(All4,bTR))*
     &           (-8*A0(MGl2) + 
     &             8*(A0(MASf2(All4,bTR)) + 
     &                B0(MB2,MGl2,MASf2(All4,bTR))*
     &                 (MGl2 - MASf2(All4,bTR))) + 
     &             MB2*(5 + 
     &                9*MGl2*
     &                 C0(MB2,0.D0,0.D0,MGl2,MASf2(All4,bTR),MGl2) + 
     &                C0(MB2,0.D0,0.D0,MGl2,MASf2(All4,bTR),
     &                  MASf2(All4,bTR))*MASf2(All4,bTR)))*
     &           UASf(All4,6,bTR))*UASfC(All4,5,bTR))/
     &      (GF*MB2**2*CKM(3,3)*CKMC(3,2)*(MGl2 - MASf2(All4,bTR)))

	ENDLOOP(All4)

#ifdef DETAILED_DEBUG
	DCONST "C8RGlu =", C8RGlu ENDL
#endif

