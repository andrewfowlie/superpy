
        cVLLSM = 1/Pi**2*
     &    (GF**2*CKM(3,3)**2*CKMC(3,2)**2*
     &      (-((D00z(MW2,MW2,Mf2(tT,3),Mf2(tT,3)) - 
     &             2*MW2*D0z(MW2,MW2,Mf2(tT,3),Mf2(tT,3)))*
     &           Mf2(tT,3)**2) + 
     &        MW2**2*(-2*(C0z(0.D0,0.D0,MW2) - 
     &              2*C0z(0.D0,MW2,Mf2(tT,3)) + 
     &              C0z(MW2,Mf2(tT,3),Mf2(tT,3)) + 
     &              12*D00z(0.D0,MW2,MW2,Mf2(tT,3)) - 
     &              6*(D00z(0.D0,0.D0,MW2,MW2) + 
     &                 D00z(MW2,MW2,Mf2(tT,3),Mf2(tT,3))) + 
     &              MW2*(D0z(0.D0,0.D0,MW2,MW2) - 
     &                 2*D0z(0.D0,MW2,MW2,Mf2(tT,3)) + 
     &                 D0z(MW2,MW2,Mf2(tT,3),Mf2(tT,3)))) + 
     &           2*(D0z(0.D0,MW2,MW2,Mf2(tT,3)) - 
     &              D0z(MW2,MW2,Mf2(tT,3),Mf2(tT,3)))*Mf2(tT,3))))

#ifdef DETAILED_DEBUG
	DCONST "cVLLSM =", cVLLSM ENDL
#endif

