(*
	EnhHbb.mod
		introduces the resummed MB into the Yukawa couplings
		this file is part of FeynHiggs
		last modified 17 May 11 th
*)


(IndexRange[Index[#]] = NoUnfold[IndexRange[Index[#]]])&/@
  {Sfermion, Chargino, Neutralino}


Conjugate[m_Mass] ^= m


ResumCoup[ c_ == rhs_ ] :=
  c == (rhs /. Mass[F[4, {g_, ___}]] -> Mdy[g]) /; FreeQ[c, S[4|6]]

ResumCoup[ other_ ] = other


M$CouplingMatrices = ResumCoup/@
  (M$CouplingMatrices /. {
    Af[t_, g1_, g2_] -> Kf[g1, g2, t]/Mass[F[t, {g1}]],
    AfC[t_, g1_, g2_] -> KfC[g1, g2, t]/Mass[F[t, {g1}]] } //.
    a_/m_Mass + b_ -> (a + m b)/m)

