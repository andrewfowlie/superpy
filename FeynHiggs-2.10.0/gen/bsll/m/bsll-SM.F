
        CALSM = 1/96.D0*
     &    (CW2*MW2*(96*D00z(0.D0,0.D0,MW2,MW2) - 
     &          96*D00z(0.D0,MW2,MW2,Mf2(tT,3)))*(MW2 - Mf2(tT,3))**2
     &         - (Mf2(tT,3)*
     &           ((72*MW2 + MB2*(12 - 16*SW2) - 12*Mf2(tT,3))*
     &              (MW2 - Mf2(tT,3)) + 
     &             4*A0(MW2)*
     &              (6*MW2 + MB2*(3 - 4*SW2) + 9*Mf2(tT,3))) - 
     &          4*MW2*A0(Mf2(tT,3))*
     &           (-MB2 + 6*MW2 + (-7 + 16*SW2)*Mf2(tT,3) + 
     &             CW2*(4*MB2 + 16*Mf2(tT,3))))/MZ2)/
     &     (CW2*(MW2 - Mf2(tT,3))**2)

#ifdef DETAILED_DEBUG
	DCONST "CALSM =", CALSM ENDL
#endif

