-(MW2*SW2*(((2*Alfa2*(1 + (1 - 2*SW2)^2))/(CW2*MW2*(-MW2 + MZ2)*SW2^2) - 
      (Alfa2*Dminus4*(1 + (1 - 2*SW2)^2))/(CW2*MW2*(-MW2 + MZ2)*SW2^2) - 
      (Alfa2*Dminus4^2*(1 + (1 - 2*SW2)^2))/(2*CW2*MW2*(-MW2 + MZ2)*SW2^2) + 
      (16*(Alfa2 - 2*Alfa2*SW2))/(CW2*MW2*(-MW2 + MZ2)*SW2^2) + 
      (8*Dminus4*(Alfa2 - 2*Alfa2*SW2))/(CW2*MW2*(-MW2 + MZ2)*SW2^2) + 
      (Dminus4^2*(Alfa2 - 2*Alfa2*SW2))/(CW2*MW2*(-MW2 + MZ2)*SW2^2))*
     A0[MW2] + ((-2*Alfa2*Dminus4^2*(MW2 + MZ2))/(CW2*MW2^2*MZ2) - 
      (2*Alfa2*MZ2*(1 + (1 - 2*SW2)^2))/(CW2*MW2^2*(-MW2 + MZ2)*SW2^2) + 
      (Alfa2*Dminus4*MZ2*(1 + (1 - 2*SW2)^2))/(CW2*MW2^2*(-MW2 + MZ2)*
        SW2^2) + (Alfa2*Dminus4^2*MZ2*(1 + (1 - 2*SW2)^2))/
       (2*CW2*MW2^2*(-MW2 + MZ2)*SW2^2) - (16*MZ2*(Alfa2 - 2*Alfa2*SW2))/
       (CW2*MW2^2*(-MW2 + MZ2)*SW2^2) - (8*Dminus4*MZ2*(Alfa2 - 2*Alfa2*SW2))/
       (CW2*MW2^2*(-MW2 + MZ2)*SW2^2) - (Dminus4^2*MZ2*(Alfa2 - 2*Alfa2*SW2))/
       (CW2*MW2^2*(-MW2 + MZ2)*SW2^2) + 
      (2*Alfa2*(MW2 + MZ2)*(10 - 20*SW2 + 4*SW2^2))/(CW2*MW2^2*MZ2*SW2^2) - 
      (Alfa2*Dminus4*(MW2 + MZ2)*(-6 + 12*SW2 + 4*SW2^2))/
       (CW2*MW2^2*MZ2*SW2^2))*A0[MZ2] + SumOver[Cha5, 2]*SumOver[Neu5, 4]*
     ((-2*Alfa2*A0[MCha2[Cha5]]*MCha2[Cha5]*VCha[Cha5, 1]*VChaC[Cha5, 1]*
        (SW*ZNeu[Neu5, 1] - CW*ZNeu[Neu5, 2])*(SW*ZNeuC[Neu5, 1] - 
         CW*ZNeuC[Neu5, 2]))/(CW2*SW2^2*(-MCha2[Cha5] + MNeu2[Neu5])*
        (-MCha2[Cha5] + MSf2[1, 1, 1])*(-MCha2[Cha5] + MSf2[1, 1, 2])) + 
      (2*Alfa2*A0[MNeu2[Neu5]]*MNeu2[Neu5]*VCha[Cha5, 1]*VChaC[Cha5, 1]*
        (SW*ZNeu[Neu5, 1] - CW*ZNeu[Neu5, 2])*(SW*ZNeuC[Neu5, 1] - 
         CW*ZNeuC[Neu5, 2]))/(CW2*SW2^2*(-MCha2[Cha5] + MNeu2[Neu5])*
        (-MNeu2[Neu5] + MSf2[1, 1, 1])*(-MNeu2[Neu5] + MSf2[1, 1, 2])) + 
      (2*Alfa2*A0[MSf2[1, 1, 1]]*MSf2[1, 1, 1]*
        ((MNeu2[Neu5] - MSf2[1, 1, 1])^(-1) + (-MCha2[Cha5] + MSf2[1, 1, 1])^
          (-1))*VCha[Cha5, 1]*VChaC[Cha5, 1]*(SW*ZNeu[Neu5, 1] - 
         CW*ZNeu[Neu5, 2])*(SW*ZNeuC[Neu5, 1] - CW*ZNeuC[Neu5, 2]))/
       (CW2*SW2^2*(-MCha2[Cha5] + MNeu2[Neu5])*(-MSf2[1, 1, 1] + 
         MSf2[1, 1, 2])) + (2*Alfa2*A0[MSf2[1, 1, 2]]*
        (MCha2[Cha5] - MNeu2[Neu5])*MSf2[1, 1, 2]*VCha[Cha5, 1]*
        VChaC[Cha5, 1]*(SW*ZNeu[Neu5, 1] - CW*ZNeu[Neu5, 2])*
        (SW*ZNeuC[Neu5, 1] - CW*ZNeuC[Neu5, 2]))/
       (CW2*SW2^2*(-MCha2[Cha5] + MNeu2[Neu5])*(MNeu2[Neu5] - MSf2[1, 1, 2])*
        (-MCha2[Cha5] + MSf2[1, 1, 2])*(-MSf2[1, 1, 1] + MSf2[1, 1, 2]))) + 
    SumOver[Cha5, 2]*SumOver[Neu5, 4]*SumOver[Sfe5, 2]*SumOver[Sfe6, 2]*
     ((-2*Alfa2*A0[MCha2[Cha5]]*MCha2[Cha5]*UCha[Cha5, 1]*UChaC[Cha5, 1]*
        USf[Sfe6, 1, 2, 1]*USfC[Sfe5, 1, 2, 2]*USfC[Sfe6, 1, 2, 1]*
        (SW*ZNeu[Neu5, 1] + CW*ZNeu[Neu5, 2])*(CB*MW*SW*USf[Sfe5, 1, 2, 2]*
          ZNeuC[Neu5, 1] + CB*CW*MW*USf[Sfe5, 1, 2, 2]*ZNeuC[Neu5, 2]))/
       (CB*CW2*MW*SW2^2*(-MCha2[Cha5] + MNeu2[Neu5])*(-MCha2[Cha5] + 
         MSf2[Sfe5, 2, 2])*(-MCha2[Cha5] + MSf2[Sfe6, 2, 1])) + 
      (2*Alfa2*A0[MNeu2[Neu5]]*MNeu2[Neu5]*UCha[Cha5, 1]*UChaC[Cha5, 1]*
        USf[Sfe6, 1, 2, 1]*USfC[Sfe5, 1, 2, 2]*USfC[Sfe6, 1, 2, 1]*
        (SW*ZNeu[Neu5, 1] + CW*ZNeu[Neu5, 2])*(CB*MW*SW*USf[Sfe5, 1, 2, 2]*
          ZNeuC[Neu5, 1] + CB*CW*MW*USf[Sfe5, 1, 2, 2]*ZNeuC[Neu5, 2]))/
       (CB*CW2*MW*SW2^2*(-MCha2[Cha5] + MNeu2[Neu5])*(-MNeu2[Neu5] + 
         MSf2[Sfe5, 2, 2])*(-MNeu2[Neu5] + MSf2[Sfe6, 2, 1])) + 
      (2*Alfa2*A0[MSf2[Sfe5, 2, 2]]*MSf2[Sfe5, 2, 2]*
        ((MNeu2[Neu5] - MSf2[Sfe5, 2, 2])^(-1) + 
         (-MCha2[Cha5] + MSf2[Sfe5, 2, 2])^(-1))*UCha[Cha5, 1]*UChaC[Cha5, 1]*
        USf[Sfe6, 1, 2, 1]*USfC[Sfe5, 1, 2, 2]*USfC[Sfe6, 1, 2, 1]*
        (SW*ZNeu[Neu5, 1] + CW*ZNeu[Neu5, 2])*(CB*MW*SW*USf[Sfe5, 1, 2, 2]*
          ZNeuC[Neu5, 1] + CB*CW*MW*USf[Sfe5, 1, 2, 2]*ZNeuC[Neu5, 2]))/
       (CB*CW2*MW*SW2^2*(-MCha2[Cha5] + MNeu2[Neu5])*(-MSf2[Sfe5, 2, 2] + 
         MSf2[Sfe6, 2, 1])) + (2*Alfa2*A0[MSf2[Sfe6, 2, 1]]*
        (MCha2[Cha5] - MNeu2[Neu5])*MSf2[Sfe6, 2, 1]*UCha[Cha5, 1]*
        UChaC[Cha5, 1]*USf[Sfe6, 1, 2, 1]*USfC[Sfe5, 1, 2, 2]*
        USfC[Sfe6, 1, 2, 1]*(SW*ZNeu[Neu5, 1] + CW*ZNeu[Neu5, 2])*
        (CB*MW*SW*USf[Sfe5, 1, 2, 2]*ZNeuC[Neu5, 1] + 
         CB*CW*MW*USf[Sfe5, 1, 2, 2]*ZNeuC[Neu5, 2]))/
       (CB*CW2*MW*SW2^2*(-MCha2[Cha5] + MNeu2[Neu5])*(MNeu2[Neu5] - 
         MSf2[Sfe6, 2, 1])*(-MCha2[Cha5] + MSf2[Sfe6, 2, 1])*
        (-MSf2[Sfe5, 2, 2] + MSf2[Sfe6, 2, 1])))))/(8*Alfa*D*Pi)
