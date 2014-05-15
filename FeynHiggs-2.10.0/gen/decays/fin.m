<< FormCalc`

lhspatt[var_, expr_] :=
  Replace[var, s_Symbol :> Pattern[s, _], {1}] -> expr

{v0, v1, col, abbr, sub} = << m/hgaZSM.amp

rc = << m/hgaZSM.rc

amp = Plus@@ v1 //. lhspatt@@@ Join[sub, abbr] //. rc

UVDivergentPart[amp] /. Re -> Identity

