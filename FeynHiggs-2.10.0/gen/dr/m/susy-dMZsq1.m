-(Alfa*CW2*MZ2)/(6*Pi*SW2) + (Alfa*(-2 + CW2/SW2 + SW2/CW2)*A0[MHp2])/
  (8*Pi) + (Alfa*(-2 + (9*CW2)/SW2 + SW2/CW2)*A0[MW2])/(8*Pi) + 
 (Alfa*(A0[MA02] + A0[Mh02] + A0[MHH2] + A0[MZ2]))/(16*CW2*Pi*SW2) + 
 (Alfa*MW2*SBA2*B0i[bb0, MZ2, Mh02, MZ2])/(4*CW2^2*Pi*SW2) + 
 (Alfa*CBA2*MW2*B0i[bb0, MZ2, MHH2, MZ2])/(4*CW2^2*Pi*SW2) + 
 ((-5*Alfa*CW2*MZ2)/(4*Pi*SW2) + (Alfa*MW2*(-(CW2/SW2) + SW2/CW2))/(2*Pi))*
  B0i[bb0, MZ2, MW2, MW2] - 
 (Alfa*SBA2*(B0i[bb00, MZ2, Mh02, MZ2] + B0i[bb00, MZ2, MHH2, MA02]))/
  (4*CW2*Pi*SW2) - (Alfa*CBA2*(B0i[bb00, MZ2, Mh02, MA02] + 
    B0i[bb00, MZ2, MHH2, MZ2]))/(4*CW2*Pi*SW2) - 
 (Alfa*(-2 + CW2/SW2 + SW2/CW2)*B0i[bb00, MZ2, MHp2, MHp2])/(4*Pi) - 
 (Alfa*(-2 + (9*CW2)/SW2 + SW2/CW2)*B0i[bb00, MZ2, MW2, MW2])/(4*Pi) - 
 (Alfa*CW2*MZ2*B0i[bb1, MZ2, MW2, MW2])/(2*Pi*SW2) + 
 ((Alfa*A0[MSf2[1, 1, Gen3]])/(8*CW2*Pi*SW2) - 
   (Alfa*B0i[bb00, MZ2, MSf2[1, 1, Gen3], MSf2[1, 1, Gen3]])/(4*CW2*Pi*SW2))*
  SumOver[Gen3, 3] + SumOver[Gen3, 3]*SumOver[Sfe3, 2]*
  ((Alfa*A0[MSf2[Sfe3, 2, Gen3]]*((1 - 2*SW2)^2*USf[Sfe3, 1, 2, Gen3]*
       USfC[Sfe3, 1, 2, Gen3] + 4*SW2^2*USf[Sfe3, 2, 2, Gen3]*
       USfC[Sfe3, 2, 2, Gen3]))/(8*CW2*Pi*SW2) + 
   (Alfa*A0[MSf2[Sfe3, 3, Gen3]]*((3 - 4*SW2)^2*USf[Sfe3, 1, 3, Gen3]*
       USfC[Sfe3, 1, 3, Gen3] + 16*SW2^2*USf[Sfe3, 2, 3, Gen3]*
       USfC[Sfe3, 2, 3, Gen3]))/(24*CW2*Pi*SW2) + 
   (Alfa*A0[MSf2[Sfe3, 4, Gen3]]*((3 - 2*SW2)^2*USf[Sfe3, 1, 4, Gen3]*
       USfC[Sfe3, 1, 4, Gen3] + 4*SW2^2*USf[Sfe3, 2, 4, Gen3]*
       USfC[Sfe3, 2, 4, Gen3]))/(24*CW2*Pi*SW2)) + 
 SumOver[Gen3, 3]*SumOver[Sfe3, 2]*SumOver[Sfe4, 2]*
  (-(Alfa*B0i[bb00, MZ2, MSf2[Sfe3, 2, Gen3], MSf2[Sfe4, 2, Gen3]]*
      ((-1 + 2*SW2)*USf[Sfe4, 1, 2, Gen3]*USfC[Sfe3, 1, 2, Gen3] + 
       2*SW2*USf[Sfe4, 2, 2, Gen3]*USfC[Sfe3, 2, 2, Gen3])*
      ((-1 + 2*SW2)*USf[Sfe3, 1, 2, Gen3]*USfC[Sfe4, 1, 2, Gen3] + 
       2*SW2*USf[Sfe3, 2, 2, Gen3]*USfC[Sfe4, 2, 2, Gen3]))/(4*CW2*Pi*SW2) - 
   (Alfa*B0i[bb00, MZ2, MSf2[Sfe3, 3, Gen3], MSf2[Sfe4, 3, Gen3]]*
     ((-3 + 4*SW2)*USf[Sfe4, 1, 3, Gen3]*USfC[Sfe3, 1, 3, Gen3] + 
      4*SW2*USf[Sfe4, 2, 3, Gen3]*USfC[Sfe3, 2, 3, Gen3])*
     ((-3 + 4*SW2)*USf[Sfe3, 1, 3, Gen3]*USfC[Sfe4, 1, 3, Gen3] + 
      4*SW2*USf[Sfe3, 2, 3, Gen3]*USfC[Sfe4, 2, 3, Gen3]))/(12*CW2*Pi*SW2) - 
   (Alfa*B0i[bb00, MZ2, MSf2[Sfe3, 4, Gen3], MSf2[Sfe4, 4, Gen3]]*
     ((-3 + 2*SW2)*USf[Sfe4, 1, 4, Gen3]*USfC[Sfe3, 1, 4, Gen3] + 
      2*SW2*USf[Sfe4, 2, 4, Gen3]*USfC[Sfe3, 2, 4, Gen3])*
     ((-3 + 2*SW2)*USf[Sfe3, 1, 4, Gen3]*USfC[Sfe4, 1, 4, Gen3] + 
      2*SW2*USf[Sfe3, 2, 4, Gen3]*USfC[Sfe4, 2, 4, Gen3]))/(12*CW2*Pi*SW2)) + 
 SumOver[Cha3, 2]*SumOver[Cha4, 2]*
  (-(Alfa*A0[MCha2[Cha4]]*((SW2*IndexDelta[Cha3, Cha4] - 
         UCha[Cha4, 1]*UChaC[Cha3, 1] - (UCha[Cha4, 2]*UChaC[Cha3, 2])/2)*
        (SW2*IndexDelta[Cha3, Cha4] - UCha[Cha3, 1]*UChaC[Cha4, 1] - 
         (UCha[Cha3, 2]*UChaC[Cha4, 2])/2) + 
       (SW2*IndexDelta[Cha3, Cha4] - VCha[Cha4, 1]*VChaC[Cha3, 1] - 
         (VCha[Cha4, 2]*VChaC[Cha3, 2])/2)*(SW2*IndexDelta[Cha3, Cha4] - 
         VCha[Cha3, 1]*VChaC[Cha4, 1] - (VCha[Cha3, 2]*VChaC[Cha4, 2])/2)))/
    (2*CW2*Pi*SW2) + (Alfa*B0i[bb00, MZ2, MCha2[Cha3], MCha2[Cha4]]*
     ((SW2*IndexDelta[Cha3, Cha4] - UCha[Cha4, 1]*UChaC[Cha3, 1] - 
        (UCha[Cha4, 2]*UChaC[Cha3, 2])/2)*(SW2*IndexDelta[Cha3, Cha4] - 
        UCha[Cha3, 1]*UChaC[Cha4, 1] - (UCha[Cha3, 2]*UChaC[Cha4, 2])/2) + 
      (SW2*IndexDelta[Cha3, Cha4] - VCha[Cha4, 1]*VChaC[Cha3, 1] - 
        (VCha[Cha4, 2]*VChaC[Cha3, 2])/2)*(SW2*IndexDelta[Cha3, Cha4] - 
        VCha[Cha3, 1]*VChaC[Cha4, 1] - (VCha[Cha3, 2]*VChaC[Cha4, 2])/2)))/
    (CW2*Pi*SW2) - (Alfa*MZ2*B0i[bb1, MZ2, MCha2[Cha3], MCha2[Cha4]]*
     ((SW2*IndexDelta[Cha3, Cha4] - UCha[Cha4, 1]*UChaC[Cha3, 1] - 
        (UCha[Cha4, 2]*UChaC[Cha3, 2])/2)*(SW2*IndexDelta[Cha3, Cha4] - 
        UCha[Cha3, 1]*UChaC[Cha4, 1] - (UCha[Cha3, 2]*UChaC[Cha4, 2])/2) + 
      (SW2*IndexDelta[Cha3, Cha4] - VCha[Cha4, 1]*VChaC[Cha3, 1] - 
        (VCha[Cha4, 2]*VChaC[Cha3, 2])/2)*(SW2*IndexDelta[Cha3, Cha4] - 
        VCha[Cha3, 1]*VChaC[Cha4, 1] - (VCha[Cha3, 2]*VChaC[Cha4, 2])/2)))/
    (2*CW2*Pi*SW2) + B0i[bb0, MZ2, MCha2[Cha3], MCha2[Cha4]]*
    ((Alfa*MCha[Cha3]*MCha[Cha4]*((SW2*IndexDelta[Cha3, Cha4] - 
          UCha[Cha4, 1]*UChaC[Cha3, 1] - (UCha[Cha4, 2]*UChaC[Cha3, 2])/2)*
         (SW2*IndexDelta[Cha3, Cha4] - VCha[Cha4, 1]*VChaC[Cha3, 1] - 
          (VCha[Cha4, 2]*VChaC[Cha3, 2])/2) + (SW2*IndexDelta[Cha3, Cha4] - 
          UCha[Cha3, 1]*UChaC[Cha4, 1] - (UCha[Cha3, 2]*UChaC[Cha4, 2])/2)*
         (SW2*IndexDelta[Cha3, Cha4] - VCha[Cha3, 1]*VChaC[Cha4, 1] - 
          (VCha[Cha3, 2]*VChaC[Cha4, 2])/2)))/(2*CW2*Pi*SW2) - 
     (Alfa*MCha2[Cha3]*((SW2*IndexDelta[Cha3, Cha4] - UCha[Cha4, 1]*
           UChaC[Cha3, 1] - (UCha[Cha4, 2]*UChaC[Cha3, 2])/2)*
         (SW2*IndexDelta[Cha3, Cha4] - UCha[Cha3, 1]*UChaC[Cha4, 1] - 
          (UCha[Cha3, 2]*UChaC[Cha4, 2])/2) + (SW2*IndexDelta[Cha3, Cha4] - 
          VCha[Cha4, 1]*VChaC[Cha3, 1] - (VCha[Cha4, 2]*VChaC[Cha3, 2])/2)*
         (SW2*IndexDelta[Cha3, Cha4] - VCha[Cha3, 1]*VChaC[Cha4, 1] - 
          (VCha[Cha3, 2]*VChaC[Cha4, 2])/2)))/(2*CW2*Pi*SW2))) + 
 SumOver[Neu3, 4]*SumOver[Neu4, 4]*
  (-(Alfa*A0[MNeu2[Neu4]]*(ZNeu[Neu4, 3]*ZNeuC[Neu3, 3] - 
       ZNeu[Neu4, 4]*ZNeuC[Neu3, 4])*(ZNeu[Neu3, 3]*ZNeuC[Neu4, 3] - 
       ZNeu[Neu3, 4]*ZNeuC[Neu4, 4]))/(8*CW2*Pi*SW2) + 
   (Alfa*B0i[bb00, MZ2, MNeu2[Neu3], MNeu2[Neu4]]*
     (ZNeu[Neu4, 3]*ZNeuC[Neu3, 3] - ZNeu[Neu4, 4]*ZNeuC[Neu3, 4])*
     (ZNeu[Neu3, 3]*ZNeuC[Neu4, 3] - ZNeu[Neu3, 4]*ZNeuC[Neu4, 4]))/
    (4*CW2*Pi*SW2) - (Alfa*MZ2*B0i[bb1, MZ2, MNeu2[Neu3], MNeu2[Neu4]]*
     (ZNeu[Neu4, 3]*ZNeuC[Neu3, 3] - ZNeu[Neu4, 4]*ZNeuC[Neu3, 4])*
     (ZNeu[Neu3, 3]*ZNeuC[Neu4, 3] - ZNeu[Neu3, 4]*ZNeuC[Neu4, 4]))/
    (8*CW2*Pi*SW2) + B0i[bb0, MZ2, MNeu2[Neu3], MNeu2[Neu4]]*
    (-(Alfa*MNeu2[Neu3]*(ZNeu[Neu4, 3]*ZNeuC[Neu3, 3] - 
         ZNeu[Neu4, 4]*ZNeuC[Neu3, 4])*(ZNeu[Neu3, 3]*ZNeuC[Neu4, 3] - 
         ZNeu[Neu3, 4]*ZNeuC[Neu4, 4]))/(8*CW2*Pi*SW2) + 
     (Alfa*MNeu[Neu3]*MNeu[Neu4]*
       (-(ZNeu[Neu4, 3]*ZNeuC[Neu3, 3] - ZNeu[Neu4, 4]*ZNeuC[Neu3, 4])^2 - 
        (ZNeu[Neu3, 3]*ZNeuC[Neu4, 3] - ZNeu[Neu3, 4]*ZNeuC[Neu4, 4])^2))/
      (16*CW2*Pi*SW2)))
