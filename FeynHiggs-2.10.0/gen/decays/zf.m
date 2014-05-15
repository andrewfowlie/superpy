<< FeynArts`

<< FormCalc`

<< FormCalc`tools`btensor`


(*SM = {"SMQCD", "EnhYuk"};*)

(*SetOptions[InsertFields, InsertionLevel -> {Particles},
(*LastSelections -> F,*)
  Restrictions -> {NoGeneration1, NoGeneration2}]*)

mrulz[t_ -> ty_] := {
  Mf[t,1] -> Mfy[ty,1],
  Mf[t,2] -> Mfy[ty,2],
  Mf[t,3] -> Mfy[ty,3],
  Mf[t,g_] -> Mfy[ty,g],
  Mf2[t,1] -> Mfy2[ty,1],
  Mf2[t,2] -> Mfy2[ty,2],
  Mf2[t,3] -> Mfy2[ty,3],
  Mf2[t,g_] -> Mfy2[ty,g],
  Mqy[t,g_] -> Mfy[ty,g]
};

Sq[Mfy[a__]] = Mfy2[a]


top = CreateTopologies[1, 1 -> 1, ExcludeTopologies -> Tadpoles]

insvirt = InsertFields[top, S[1] -> S[1], Model -> SM]

virt = CreateFeynAmp[insvirt]

top = CreateCTTopologies[1, 1 -> 1, ExcludeTopologies -> TadpoleCTs]

insct = InsertFields[top, S[1] -> S[1], Model -> SM]

ct = CreateFeynAmp[insct]


seH = Plus@@ CalcFeynAmp[virt, ct, OnShell -> False] //.
  Abbr[] /. Pair[_k, _k] -> K2 /.
  Finite -> 1 /.
  _Enh -> 1

rcs = CalcRenConst[seH]

seHR = seH /. rcs

dseHR = D[seHR, K2] /. K2 -> MH2 (*/. Re -> Identity*)

final[x_] := SplitSums[x /. mrulz[3 -> tD] /. mrulz[4 -> bD] /.
  ToOldBRules /.
  { _Enh -> 1,
    Finite -> 1,
    Sqrt[2] -> sqrt2,
    1/Sqrt[2] -> 1/sqrt2 }]



Simplify[final[dseHR] //.
  a_ x_ + b_ Re[x_] :> a I Im[x] /; a + b == 0] //.
  a_ x_ + b_ x_ :> x (a + b)

