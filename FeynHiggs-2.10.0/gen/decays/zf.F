
        dse_hh = Alfa/(Pi*SW2)*
     &   (1/CW2*(-(cI/16.D0*Im(B0(MH2,MZ2,MZ2))) + 
     &        cI/8.D0*(Im(B1(MH2,MZ2,MZ2)) + 
     &           MH2*Im(DB1(MH2,MZ2,MZ2)))) - 
     &     cI/8.D0*Im(B0(MH2,MW2,MW2)) + 
     &     (9*cI)/(32.D0*MW2)*(MH2**2*Im(DB0(MH2,MH2,MH2))) + 
     &     cI*((-(1/8.D0*MH2) + 1/(16.D0*MW2)*(MH2**2 + 12*MW2**2))*
     &         Im(DB0(MH2,MW2,MW2)) + 
     &        (-(1/(16.D0*CW2)*MH2) - 
     &           1/(32.D0*CW2**2*MW2)*
     &            (-16*MW2**2 + CW2*(4*MW2*MZ2 - CW2*MH2**2)))*
     &         Im(DB0(MH2,MZ2,MZ2))) + 
     &     cI/4.D0*(Im(B1(MH2,MW2,MW2)) + MH2*Im(DB1(MH2,MW2,MW2))))

	LOOP(Gen3, 1,3,1)

        dse_hh = dse_hh + 
     &   Alfa/(MW2*Pi*SW2)*
     &    (-(cI/4.D0*((Im(B1(MH2,Mf2(2,Gen3),Mf2(2,Gen3))) + 
     &             MH2*Im(DB1(MH2,Mf2(2,Gen3),Mf2(2,Gen3))))*
     &           Mf2(2,Gen3))) - 
     &      cI/2.D0*(Im(DB0(MH2,Mf2(2,Gen3),Mf2(2,Gen3)))*
     &         Mf2(2,Gen3)**2) - 
     &      (3*cI)/4.D0*((Im(B1(MH2,Mf2(3,Gen3),Mf2(3,Gen3))) + 
     &            MH2*Im(DB1(MH2,Mf2(3,Gen3),Mf2(3,Gen3))))*
     &          Mf2(3,Gen3) + 
     &         (Im(B1(MH2,Mf2(4,Gen3),Mf2(4,Gen3))) + 
     &            MH2*Im(DB1(MH2,Mf2(4,Gen3),Mf2(4,Gen3))))*
     &          Mf2(4,Gen3)) - 
     &      (3*cI)/2.D0*(Im(DB0(MH2,Mf2(3,Gen3),Mf2(3,Gen3)))*
     &          Mf2(3,Gen3)**2 + 
     &         Im(DB0(MH2,Mf2(4,Gen3),Mf2(4,Gen3)))*Mf2(4,Gen3)**2)
     &      )

	ENDLOOP(Gen3)

