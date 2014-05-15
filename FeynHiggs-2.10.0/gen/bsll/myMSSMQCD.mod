LoadModel["FVMSSM"]
(*LoadModel["MSSMQCD"]*)


(IndexRange[Index[#]] = NoUnfold[IndexRange[Index[#]]])&/@
  {Sfermion, Chargino, Neutralino}


ResumCoup[ c_ == rhs_ ] :=
  c == (rhs /. Mass[F[4, {g_, ___}]] -> Mdy[g]) /; FreeQ[c, S[4|6]]

ResumCoup[ other_ ] = other


EnhCoup[ (c:C[F[4, {j1_, _}], -F[4, _], S[h:1|2]]) == rhs_ ] :=
  c == Hff[h, j1] rhs

EnhCoup[ other_ ] = other


M$CouplingMatrices = ResumCoup/@ EnhCoup/@ M$CouplingMatrices

