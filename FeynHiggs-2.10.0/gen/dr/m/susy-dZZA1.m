(2*(-(Alfa*(-CW2 + SW2)*A0[MHp2])/(4*CW*Pi*SW) - 
   (Alfa*(-4*CW - CW2/CW + SW2/CW)*A0[MW2])/(4*Pi*SW) - 
   (Alfa*MW2*(CW/SW + SW/CW)*B0i[bb0, 0, MW2, MW2])/(2*Pi) + 
   (Alfa*(-CW2 + SW2)*B0i[bb00, 0, MHp2, MHp2])/(2*CW*Pi*SW) + 
   (Alfa*(-4*CW - CW2/CW + SW2/CW)*B0i[bb00, 0, MW2, MW2])/(2*Pi*SW) + 
   SumOver[Gen3, 3]*SumOver[Sfe3, 2]*
    (-(Alfa*A0[MSf2[Sfe3, 2, Gen3]]*((-1 + 2*SW2)*USf[Sfe3, 1, 2, Gen3]*
          USfC[Sfe3, 1, 2, Gen3] + 2*SW2*USf[Sfe3, 2, 2, Gen3]*
          USfC[Sfe3, 2, 2, Gen3]))/(4*CW*Pi*SW) + 
     (Alfa*B0i[bb00, 0, MSf2[Sfe3, 2, Gen3], MSf2[Sfe3, 2, Gen3]]*
       ((-1 + 2*SW2)*USf[Sfe3, 1, 2, Gen3]*USfC[Sfe3, 1, 2, Gen3] + 
        2*SW2*USf[Sfe3, 2, 2, Gen3]*USfC[Sfe3, 2, 2, Gen3]))/(2*CW*Pi*SW) - 
     (Alfa*A0[MSf2[Sfe3, 3, Gen3]]*((-3 + 4*SW2)*USf[Sfe3, 1, 3, Gen3]*
         USfC[Sfe3, 1, 3, Gen3] + 4*SW2*USf[Sfe3, 2, 3, Gen3]*
         USfC[Sfe3, 2, 3, Gen3]))/(6*CW*Pi*SW) + 
     (Alfa*B0i[bb00, 0, MSf2[Sfe3, 3, Gen3], MSf2[Sfe3, 3, Gen3]]*
       ((-3 + 4*SW2)*USf[Sfe3, 1, 3, Gen3]*USfC[Sfe3, 1, 3, Gen3] + 
        4*SW2*USf[Sfe3, 2, 3, Gen3]*USfC[Sfe3, 2, 3, Gen3]))/(3*CW*Pi*SW) - 
     (Alfa*A0[MSf2[Sfe3, 4, Gen3]]*((-3 + 2*SW2)*USf[Sfe3, 1, 4, Gen3]*
         USfC[Sfe3, 1, 4, Gen3] + 2*SW2*USf[Sfe3, 2, 4, Gen3]*
         USfC[Sfe3, 2, 4, Gen3]))/(12*CW*Pi*SW) + 
     (Alfa*B0i[bb00, 0, MSf2[Sfe3, 4, Gen3], MSf2[Sfe3, 4, Gen3]]*
       ((-3 + 2*SW2)*USf[Sfe3, 1, 4, Gen3]*USfC[Sfe3, 1, 4, Gen3] + 
        2*SW2*USf[Sfe3, 2, 4, Gen3]*USfC[Sfe3, 2, 4, Gen3]))/(6*CW*Pi*SW)) + 
   SumOver[Cha3, 2]*((Alfa*A0[MCha2[Cha3]]*(4*SW2 - 2*UCha[Cha3, 1]*
         UChaC[Cha3, 1] - UCha[Cha3, 2]*UChaC[Cha3, 2] - 
        2*VCha[Cha3, 1]*VChaC[Cha3, 1] - VCha[Cha3, 2]*VChaC[Cha3, 2]))/
      (4*CW*Pi*SW) - (Alfa*B0i[bb00, 0, MCha2[Cha3], MCha2[Cha3]]*
       (4*SW2 - 2*UCha[Cha3, 1]*UChaC[Cha3, 1] - UCha[Cha3, 2]*
         UChaC[Cha3, 2] - 2*VCha[Cha3, 1]*VChaC[Cha3, 1] - 
        VCha[Cha3, 2]*VChaC[Cha3, 2]))/(2*CW*Pi*SW))))/MZ2
