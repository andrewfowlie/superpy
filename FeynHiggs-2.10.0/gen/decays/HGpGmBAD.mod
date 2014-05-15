(*
	HGpGm.mod
	Modification of the Higgs -- G+ -- G- coupling
	for the H->gamma gamma decay
*)

MHCoup[ (c:C[S[1], _. S[6], _. -S[6]]) == _ ] :=
  c == {-I EL/(2 MW SW) SBA MHiggs2[1]}

MHCoup[ (c:C[S[2], _. S[6], _. -S[6]]) == _ ] :=
  c == {-I EL/(2 MW SW) CBA MHiggs2[2]}

MHCoup[other_] := other


M$CouplingMatrices = MHCoup/@ M$CouplingMatrices

