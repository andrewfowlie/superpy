
	C8LSM = 0

	LOOP(Gen4, 1,3,1)

        dup1 = 2*(MB2 - MW2)*MW2 + 
     &    Mf2(3,Gen4)*(-MB2 + MW2 + Mf2(3,Gen4))

        C8LSM = C8LSM + 
     &    1/4.D0*(CKM(Gen4,3)*CKMC(Gen4,2)*
     &        (dup1*(-2*A0(MW2) + 2*A0(Mf2(3,Gen4))) - 
     &          (MW2 - Mf2(3,Gen4))*
     &           (-2*(MB2*MW2 + dup1*B0(MB2,MW2,Mf2(3,Gen4))) + 
     &             MB2*(-1 + 
     &                2*C0(MB2,0.D0,0.D0,MW2,Mf2(3,Gen4),Mf2(3,Gen4))*
     &                 (MB2 - 2*MW2 - Mf2(3,Gen4)))*Mf2(3,Gen4))))/
     &      (MB2**2*CKM(3,3)*CKMC(3,2)*(MW2 - Mf2(3,Gen4)))

	ENDLOOP(Gen4)

#ifdef DETAILED_DEBUG
	DCONST "C8LSM =", C8LSM ENDL
#endif

