(*
	ResumHbb.mod
		adds variables for KW's resummation factors
		to the Hbb coupling
		this file is part of FeynHiggs
		last modified 23 Jun 11 th
*)


EnhCoup[ (c:C[_. F[4, {j1_, _}],
              _. F[4, {j2_, _}],
              _. S[h:1|2|3]]) == {cL_, cR_} ] :=
  c == {HffDb[sub1L,h,4,j1] cL, Conjugate[HffDb[sub1L,h,4,j2]] cR}

EnhCoup[ (c:C[F[15, _], _. F[4, _], _. S[14, _]]) == rhs_ ] :=
  c == (rhs /. GS -> GSDb)

EnhCoup[other_] = other


M$CouplingMatrices = EnhCoup/@ M$CouplingMatrices

