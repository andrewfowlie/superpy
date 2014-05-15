(*
	HMixExt.mod
		Add-on model file which adds two new particles:
		S[0, {h}] = Sum[UHiggs[h,i] S[i], {i, 3}]
		S[10, {h}] = Sum[ZHiggs[h,i] S[i], {i, 3}]
		this file is part of FeynArts
		last modified 17 Mar 09 th
*)


IndexRange[ Index[Higgs] ] = Range[3]

M$ClassesDescription = Flatten[{
  M$ClassesDescription,
  S[0] == {
	SelfConjugate -> True,
	Indices -> {Index[Higgs]},
	InsertOnly -> {External},
	Mass -> MHiggs,
	PropagatorLabel -> ComposedChar["H", Index[Higgs]],
	PropagatorType -> ScalarDash,
	PropagatorArrow -> None },
  S[10] == {
	SelfConjugate -> True,
	Indices -> {Index[Higgs]},
	InsertOnly -> {External},
	Mass -> MHiggs,
	PropagatorLabel -> ComposedChar["H", Index[Higgs], Null, "\\hat"],
	PropagatorType -> ScalarDash,
	PropagatorArrow -> None }
}]


Block[ {selH, sign, addZ, NewCoup, UZPerm, coup, oldcoup, newcoup},

selH[S[h:1|2|3]] := {S[h], S[10, h]};
selH[f_] := {f};

addZ[c_C] := Block[ {hi = 0},
  add[Flatten/@ Transpose[addZ/@ Sort[List@@ c]], sign[c]] ];
addZ[S[10, h_]] :=
  {S[10, {#}], #, h}& @ ToExpression["h" <> ToString[++hi]];
addZ[f_] := {f, {}, {}};

NewCoup[c_ == rhs_] := {c == rhs} /; FreeQ[c, S[1|2|3]];

NewCoup[c_ == rhs_] := {c == rhs,
Block[ {sign, add, comb, coup, AddTo},
  sign = If[ToGeneric[c] === C[S, S, V], Signature, 1 &];
  add[{{f__}, i_, j_}, s_] := coup[f] += s rhs *
    (Plus@@ (sign[#] Times@@ MapThread[ZHiggs, {i, #}]&)/@
      Permutations[j]);
  comb = Flatten[Outer[C, Sequence@@ selH/@ c]];
  Union[ addZ/@ Select[comb, !FreeQ[#, S[10, _]]&],
    SameTest -> (#1[[1]] === #2[[1]] &) ]
]};


UZPerm[_[_[c__]], rhs_] :=
  {#, # /. {S[10, h_] -> S[0, h], ZHiggs -> UHiggs}}& @
    (C[c] == (*Simplify @*) rhs);


_coup = 0;
oldcoup = First/@ NewCoup/@ M$CouplingMatrices;
_coup =.;
newcoup = Apply[UZPerm, DownValues[coup], 1];


If[ TrueQ[$JustNewCouplings], oldcoup = {} ];

M$CouplingMatrices = Flatten[{oldcoup, newcoup}];

]

