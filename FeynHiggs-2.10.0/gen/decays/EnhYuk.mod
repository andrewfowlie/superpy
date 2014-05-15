(*
	EnhYuk.mod
		adds enhancement factors and introduces the
		resummed MB into the Yukawa couplings
		this file is part of FeynHiggs
		last modified 18 May 11 th
*)


Conjugate[m_Mass] ^= m


ResumCoup[ c_ == rhs_ ] :=
  c == (rhs /. Mass[F[t:3|4, {g_, ___}]] -> Mqy[t, g]) /; FreeQ[c, S[4|6]]

ResumCoup[ other_ ] = other


Enh/: Enh[h__, {t1_, j1_}, {t2_, j2_}] IndexDelta[j1_, j2_] :=
  Enh[h, {t1, j1}, {t2, j1}] IndexDelta[j1, j2]

EnhCoup[ (c:C[_. F[t1:3|4, {j1_, _}],
              _. F[t2:3|4, {j2_, _}],
              _. S[h_]]) == rhs_ ] :=
  c == Enh[Q, h, {t1, j1}, {t2, j2}] rhs

EnhCoup[ (c:C[_. S[h:1|2|3|4|5|6],
              _. S[t1:13|14, {___, j1_, _}],
              _. S[t2:13|14, {___, j2_, _}]]) == rhs_ ] :=
  c == Enh[Sq, h, {t1, j1}, {t2, j2}] rhs

EnhCoup[ (c:C[_. S[h1:1|2|3|4|5|6],
              _. S[h2:1|2|3|4|5|6],
              _. S[t1:13|14, {___, j1_, _}],
              _. S[t2:13|14, {___, j2_, _}]]) == rhs_ ] :=
  c == Enh[Sq, h1, h2, {t1, j1}, {t2, j2}] rhs

EnhCoup[other_] = other


M$CouplingMatrices = EnhCoup/@ ResumCoup/@ M$CouplingMatrices


Scan[
  (DownValues[#] = DeleteCases[DownValues[#], _[_, _Symbol]])&,
  {MLE, MQU, MQD, Mf, Mf2} ]

