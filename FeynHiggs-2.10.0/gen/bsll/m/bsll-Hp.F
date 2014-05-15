
        CALHp = -(1/24.D0*
     &      (((-1 + 4*CW2)*MB*TB2*Mf(bTR,3) + 3*Mf2(tT,3))*
     &         (MHp2*(A0(Mf2(tT,3)) - Mf2(tT,3)) + 
     &           Mf2(tT,3)*(-A0(MHp2) + Mf2(tT,3))))/
     &       (CW2*MZ2*TB2*(MHp2 - Mf2(tT,3))**2))

#ifdef DETAILED_DEBUG
	DCONST "CALHp =", CALHp ENDL
#endif

        CARHp = -(1/24.D0*
     &      (Mf(bTR,2)*(4*MB*SW2 + 3*TB2*Mf(bTR,3))*
     &         (-(MHp2*A0(Mf2(tT,3))) + 
     &           (MHp2 + A0(MHp2) - Mf2(tT,3))*Mf2(tT,3)))/
     &       (CW2*MZ2*(MHp2 - Mf2(tT,3))**2))

#ifdef DETAILED_DEBUG
	DCONST "CARHp =", CARHp ENDL
#endif

