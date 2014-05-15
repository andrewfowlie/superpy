(*
	KW.mod
		model file which adds *only* a few counter terms to MSSM
		by Karina Williams
		last modified 10 Jun 11 th
*)


Keeptstbsb = ExcludeParticles -> {
  _U, _V,
  S[0 | 1 | 2 | 3 | 4 | 5 | 6 | 11 | 12], 
  S[13, {_, 1 | 2, _}], S[14, {_, 1 | 2, _}],
  F[1 | 2 | 11 | 12], F[3, {1 | 2, _}], F[4, {1 | 2, _}] }
	
KeeptstbsbGZ = ExcludeParticles -> {
  _U, V[1 | 3],
  S[0 | 1 | 2 | 3 | 5 | 6 | 11 | 12 ], 
  S[13, {_, 1 | 2, _}], S[14, {_, 1 | 2, _}],
  F[1 | 2 | 11 | 12], F[3, {1 | 2, _}], F[4, {1 | 2, _}] }


epGP = SqrtEGl;
emGP = Conjugate[SqrtEGl]


newC/: (newC[lhs__] == rhs_) :=
Block[ {pos = Position[M$CouplingMatrices, C[lhs]]},
  If[ Length[pos] =!= 0,
    M$CouplingMatrices[[ pos[[1, 1]], 2 ]] = rhs /.
      oldC[r___] :> (M$CouplingMatrices[[ pos[[1, 1]], 2 ]] /. {r}),
    AppendTo[M$CouplingMatrices, C[lhs] == rhs]
  ]
]


newC[S[1], S[1]] == -I {
  {0, dZh0h01},
  {0, dMh0h0sq1 + Mh0tree2 dZh0h01}
}

newC[S[1], S[2]] == -I {
  {0, dZh0HH1},
  {0, dMh0HHsq1 + (Mh0tree2 + MHHtree2)/2 dZh0HH1}
}

newC[S[1], S[3]] == -I {
  {0, dZh0A01},
  {0, dMh0A0sq1 + (Mh0tree2 + MA0tree2)/2 dZh0A01}
}

newC[S[2], S[2]] == -I {
  {0, dZHHHH1},
  {0, dMHHHHsq1 + MHHtree2 dZHHHH1}
}

newC[S[2], S[3]] == -I {
  {0, dZHHA01},
  {0, dMHHA0sq1 + (MHHtree2 + MA0tree2)/2 dZHHA01}
}

newC[S[3], S[3]] == -I {
  {0, dZA0A01},
  {0, dMA0A0sq1 + MA0tree2 dZA0A01}
}

newC[S[1], S[4]] == -I {
  {0, 0},
  {0, dMh0G0sq1}
}

newC[S[2], S[4]] == -I {
  {0, 0},
  {0, dMHHG0sq1}
}

newC[S[3], S[4]] == -I {
  {0, dZA0G01},
  {0, dMA0G0sq1 + MA0tree2 dZA0G01/2}
}

newC[S[3], V[2]] == 1/2 {
  {0, -MZ dZA0G01},
  {0,  MZ dZA0G01}
}


sss[f__] := Cases[M$CouplingMatrices, (# == x_) -> x][[1,1,1]]&[ S/@ C[f] ];
ch0h0h0 = sss[1, 1, 1];
ch0h0HH = sss[1, 1, 2];
ch0HHHH = sss[1, 2, 2];
cHHHHHH = sss[2, 2, 2];
ch0A0A0 = sss[1, 3, 3];
cHHA0A0 = sss[2, 3, 3];
ch0A0G0 = sss[1, 3, 4];
cHHA0G0 = sss[2, 3, 4]

dEW1 = dSW1/SW (SW^2 - CW^2)/CW^2 + dMZsq1/(2 MZ^2) + dZe1

newC[S[1], S[1], S[1]] == ch0h0h0 {{1,
  3/2 dZh0h01 +
    3/2 dZh0HH1 ch0h0HH/ch0h0h0 +
    dEW1 + SB CB dTB1 cHHHHHH/ch0h0h0
}}

newC[S[2], S[2], S[2]] == cHHHHHH {{1,
  3/2 dZHHHH1 +
    3/2 dZh0HH1 ch0HHHH/cHHHHHH +
    dEW1 - SB CB dTB1 ch0h0h0/cHHHHHH
}}

newC[S[1], S[1], S[2]] == ch0h0HH {{1,
  dZh0h01 + dZHHHH1/2 +
    dZh0HH1/2 ch0h0h0/ch0h0HH + dZh0HH1 ch0HHHH/ch0h0HH +
    dEW1 - SB CB dTB1 ch0HHHH/ch0h0HH
}}

newC[S[1], S[2], S[2]] == ch0HHHH {{1,
  dZHHHH1 + dZh0h01/2 +
    dZh0HH1/2 cHHHHHH/ch0HHHH + dZh0HH1 ch0h0HH/ch0HHHH +
    dEW1 + SB CB dTB1 ch0h0HH/ch0HHHH
}}

newC[S[1], S[3], S[3]] == ch0A0A0 {{1,
  dZA0A01 + dZh0h01/2 +
    dZh0HH1/2 cHHA0A0/ch0A0A0 + dZA0G01 ch0A0G0/ch0A0A0 +
    dEW1 - SB CB dTB1 cHHA0A0/ch0A0A0
}}

newC[S[2], S[3], S[3]] == cHHA0A0 {{1,
  dZHHHH1/2 + dZA0A01 +
    dZh0HH1/2 ch0A0A0/cHHA0A0 + dZA0G01 cHHA0G0/cHHA0A0 +
    dEW1 + SB CB dTB1 ch0A0A0/cHHA0A0
}}


(* added 29 Nov 07 th *)

newC[ S[1], V[2], V[1] ] == I EL MW/(SW CW^2) {
  {0, dZZA1/2}
}


ctFFS[1, 2] := dZh0h01/2 - dZh0HH1/2 CA/SA + dTB1 SB^2;
ctFFS[2, 2] := dZHHHH1/2 - dZh0HH1/2 SA/CA + dTB1 SB^2;
ctFFS[3, 2] := dZA0A01/2 - dZA0G01/2 1/TB  + dTB1 SB^2

ctFFS[1, 3] := dZh0h01/2 + dZh0HH1/2 SA/CA - dTB1 CB^2;
ctFFS[2, 3] := dZHHHH1/2 + dZh0HH1/2 CA/SA - dTB1 CB^2;
ctFFS[3, 3] := dZA0A01/2 + dZA0G01/2 TB    - dTB1 CB^2

ctFFS[1, 4] := dZh0h01/2 - dZh0HH1/2 CA/SA + dTB1 SB^2;
ctFFS[2, 4] := dZHHHH1/2 - dZh0HH1/2 SA/CA + dTB1 SB^2;
ctFFS[3, 4] := dZA0A01/2 - dZA0G01/2 1/TB  + dTB1 SB^2

ctFFSLR[h_, t_, j1_, j2_] := IndexDelta[j1, j2] Mass[F[t, {j1}]] *
  {1, dZe1 - dMWsq1/(2 MW^2) - dSW1/SW +
        dMf1[t, j1]/Mass[F[t, {j1}]] +
        ctFFS[h, t]}

mdZfLR1[ type_, j1_, j2_ ] :=
  Mass[F[type, {j1}]]/2 dZfL1[type, j1, j2] +
    Mass[F[type, {j2}]]/2 Conjugate[dZfR1[type, j2, j1]]

mdZfRL1[ type_, j1_, j2_ ] :=
  Mass[F[type, {j1}]]/2 dZfR1[type, j1, j2] +
    Mass[F[type, {j2}]]/2 Conjugate[dZfL1[type, j2, j1]]

ctFFSL[h_, t_, j1_, j2_] :=
  ctFFSLR[h, t, j1, j2] + {0, mdZfLR1[t, j1, j2]}

ctFFSR[h_, t_, j1_, j2_] :=
  ctFFSLR[h, t, j1, j2] + {0, mdZfRL1[t, j1, j2]}


newC[F[2, {j1}], -F[2, {j2}], S[1]] == I/2 EL SA/(CB MW SW) {
  ctFFSL[1, 2, j1, j2],
  ctFFSR[1, 2, j1, j2]
}

newC[F[2, {j1}], -F[2, {j2}], S[2]] == -I/2 EL CA/(CB MW SW) {
  ctFFSL[2, 2, j1, j2],
  ctFFSR[2, 2, j1, j2]
}

newC[F[2, {j1}], -F[2, {j2}], S[3]] == EL TB/(2 MW SW) {
  ctFFSL[3, 2, j1, j2],
 -ctFFSR[3, 2, j1, j2]
}

newC[F[3, {j1, o1}], -F[3, {j2, o2}], S[1]] == -I/2 EL CA/(MW SB SW) *
    IndexDelta[o1, o2] {
  ctFFSL[1, 3, j1, j2],
  ctFFSR[1, 3, j1, j2]
}

newC[F[3, {j1, o1}], -F[3, {j2, o2}], S[2]] == -I/2 EL SA/(MW SB SW) *
    IndexDelta[o1, o2] {
  ctFFSL[2, 3, j1, j2],
  ctFFSR[2, 3, j1, j2]
}

newC[F[3, {j1, o1}], -F[3, {j2, o2}], S[3]] == EL/(2 MW SW TB) *
    IndexDelta[o1, o2] {
  ctFFSL[3, 3, j1, j2],
 -ctFFSR[3, 3, j1, j2]
}

newC[F[4, {j1, o1}], -F[4, {j2, o2}], S[1]] == I/2 EL SA/(CB MW SW) *
    IndexDelta[o1, o2] {
  ctFFSL[1, 4, j1, j2],
  ctFFSR[1, 4, j1, j2]
}

newC[F[4, {j1, o1}], -F[4, {j2, o2}], S[2]] == -I/2 EL CA/(CB MW SW) *
    IndexDelta[o1, o2] {
  ctFFSL[2, 4, j1, j2],
  ctFFSR[2, 4, j1, j2]
}

newC[F[4, {j1, o1}], -F[4, {j2, o2}], S[3]] == EL TB/(2 MW SW) *
    IndexDelta[o1, o2] {
  ctFFSL[3, 4, j1, j2],
 -ctFFSR[3, 4, j1, j2]
}

newC[F[15, {g1}], -F[3, {j1, o1}], S[13, {s2, j2, o2}]] == I GS *
    Sqrt[2] SUNT[g1, o1, o2] IndexDelta[j1, j2] {
  { emGP Conjugate[USf[3, j1][s2, 2]]},
  {-epGP Conjugate[USf[3, j1][s2, 1]]}
}

newC[F[15, {g1}], -F[4, {j1, o1}], S[14, {s2, j2, o2}]] == I GS *
    Sqrt[2] SUNT[g1, o1, o2] IndexDelta[j1, j2] {
  { emGP Conjugate[USf[4, j1][s2, 2]]},
  {-epGP Conjugate[USf[4, j1][s2, 1]]}
}

newC[F[15, {g1}], F[3, {j1, o1}], -S[13, {s2, j2, o2}]] == I GS *
    Sqrt[2] SUNT[g1, o2, o1] IndexDelta[j1, j2] {
  {-emGP USf[3, j1][s2, 1]},
  { epGP USf[3, j1][s2, 2]}
}

newC[F[15, {g1}], F[4, {j1, o1}], -S[14, {s2, j2, o2}]] == I GS *
    Sqrt[2] SUNT[g1, o2, o1] IndexDelta[j1, j2] {
  {-emGP USf[4, j1][s2, 1]},
  { epGP USf[4, j1][s2, 2]}
}

newC[-F[4, {j1, o1}], F[4, {j2, o2}]] == I IndexDelta[o1, o2] {
  {0, -AddHC[dZfL1[4, j1, j2]]},
  {0,  AddHC[dZfR1[4, j1, j2]]},
  {0, -mdZfLR1[4, j1, j2] - IndexDelta[j1, j2] dMf1[4, j1]},
  {0, -mdZfRL1[4, j1, j2] - IndexDelta[j1, j2] dMf1[4, j1]}
}

(***************************)      
            
C2BA = CBA^2 - SBA^2;
S2BA = 2 SBA CBA;
C2AB = CAB^2 - SAB^2;
S2AB = 2 SAB CAB;

RenConst[ dMh0h0sq1 ] := dMA0A0sq1 CBA^2 + dMZsq1 SAB^2 + 
  EL (dTHH1 CBA SBA^2 - dTh01 SBA (1 + CBA^2))/(2 MZ SW CW) -
  dTB1 S2B (MA0tree2 S2BA - MZ2 S2AB)/2;
RenConst[ dMh0HHsq1 ] := -dMA0A0sq1 S2BA/2 - dMZsq1 S2AB/2 + 
  EL (-dTHH1 SBA^3 - dTh01 CBA^3)/(2 MZ SW CW) -
  dTB1 S2B (MA0tree2 C2BA + MZ^2 C2AB)/2;
RenConst[ dMHHHHsq1 ] := dMA0A0sq1 SBA^2 + dMZsq1 CAB^2 -
  EL (dTHH1 CBA (1 + SBA^2) - dTh01 SBA CBA^2)/(2 MZ SW CW) +
  dTB1 S2B (MA0tree2 S2BA - MZ^2 S2AB)/2

RenConst[ dMh0A0sq1 ] := EL (-dTA01 SBA)/(2 MZ SW CW);
RenConst[ dMHHA0sq1 ] := EL (-dTA01 CBA)/(2 MZ SW CW);
RenConst[ dMA0A0sq1 ] := dMHpsq1 - dMWsq1

RenConst[ dMh0G0sq1 ] := -dMHHA0sq1;
RenConst[ dMHHG0sq1 ] := dMh0A0sq1;
RenConst[ dMA0G0sq1 ] := EL/(2 MZ SW CW) (dTHH1 SBA - dTh01 CBA) -
  dTB1 SB CB MA0tree2

RenConst[ dTh01 ] := TadpoleRC[S[1]];
RenConst[ dTHH1 ] := TadpoleRC[S[2]];
RenConst[ dTA01 ] := TadpoleRC[S[3]]

RenConst[ dZh0h01 ] := UVDivergentPart[FieldRC[S[1]]];
RenConst[ dZHHHH1 ] := UVDivergentPart[FieldRC[S[2]]];
RenConst[ dZA0A01 ] := UVDivergentPart[FieldRC[S[3]]];
RenConst[ dZG0G01 ] := UVDivergentPart[FieldRC[S[4]]]

RenConst[ dZh0HH1 ] := SA CA (dZh0h01 - dZHHHH1)/C2A;
RenConst[ dZh0A01 ] := 0;
RenConst[ dZHHA01 ] := 0

RenConst[ dZh0G01 ] := 0;
RenConst[ dZHHG01 ] := 0;
RenConst[ dZA0G01 ] := SB CB (dZh0h01 - dZHHHH1)/C2A

RenConst[ dMWsq1 ] := MassRC[V[3]];
RenConst[ dMZsq1 ] := MassRC[V[2]];
RenConst[ dMHpsq1 ] := MassRC[S[5]]

RenConst[ dTB1 ] := (dZh0h01 - dZHHHH1)/(2 C2A);
RenConst[ dSW1 ] := CW2/SW/2 (dMZsq1/MZ^2 - dMWsq1/MW^2)

(*RenConst[ dEL1 ] := EL (dSW1/SW - 
  (ReTilde[SelfEnergy[V[3] -> V[3], 0]] - dMWsq1)/(2 MW^2))*)
(*RenConst[ dEL1 ] := EL dZe1*)


(* nb for gluino contrib to bottom mass, *)
(*  both LVectorCoeff[sff] and  RVectorCoeff[sff] *)
(* are real *)

RenConst[ dZfL1[t_, j1_, j1_] ] :=
Block[ {m = TheMass[F[t, {j1}]], sff = SelfEnergy[F[t, {j1}]]},
  FieldRC[F[t, {j1}]][[1]] +
    1/(2 m) ReTilde[LScalarCoeff[sff] - RScalarCoeff[sff]]
]

RenConst[ dZfR1[t_, j1_, j1_] ] :=
Block[ {m = TheMass[F[t, {j1}]], sff = SelfEnergy[F[t, {j1}]]},
  FieldRC[F[t, {j1}]][[2]] -
    1/(2 m) ReTilde[LScalarCoeff[sff] - RScalarCoeff[sff]]
]

