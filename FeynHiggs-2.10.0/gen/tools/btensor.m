(*
	btensor.m
		explicit decompositions of the two-point
		tensor-coefficient functions
		this file is part of FormCalc
		last modified 10 Nov 09 th
*)


ok[0] = True;
ok[MW2] = True;
ok[MZ2] = True;
ok[MT2] = True;
ok[MH2] = True;
ok[Mh02] = True;
ok[MHH2] = True;
ok[MA02] = True;
ok[MHp2] = True;
(*
ok[_MSf2] = True;
ok[_MNeu2] = True;
ok[_MCha2] = True;
*)
ok[_] = False;
ok[m__] := VectorQ[{m}, ok]

A0[0] = 0


B0[0, 0, 0] = BAD[B0]	(* divergent, must cancel *)

B0[0, m_, m_] := A0[m]/m - 1 /; m =!= 0

B0[m_, 0, m_] := A0[m]/m + 1 /; m =!= 0

B0[m_, m_, 0] := A0[m]/m + 1 /; m =!= 0

B0[0, m1_, m2_] := (A0[m1] - A0[m2])/(m1 - m2) /; ok[m1, m2]


DeltaB0[p_, m1_, m2_] := (B0[p, m1, m2] - B0[0, m1, m2])/p

DeltaB1[p_, m1_, m2_] := (B1[p, m1, m2] - B1[0, m1, m2])/p


B1[0, m_, m_] := -B0[0, m, m]/2

B1[0, m1_, m2_] := -B0[0, m1, m2]/2 - (m1 - m2)/2 DB0[0, m1, m2]

B1[p_, m1_, m2_] := -B0[p, m1, m2]/2 - (m1 - m2)/2 DeltaB0[p, m1, m2]


B00[0, m1_, m2_] := 1/4 (A0[m2] + m1 B0[0, m1, m2] + (m1 + m2)/2)

B00[p_, m1_, m2_] :=
  -(p - 3 (m1 + m2))/18 + m1 B0[p, m1, m2]/3 +
  (A0[m2] + (p + m1 - m2) B1[p, m1, m2])/6


B11[0, m1_, m2_] := 1/18 - 2/3 B1[0, m1, m2] -
  2/3 (m1 - m2) DB1[0, m1, m2] - 1/3 m1 DB0[0, m1, m2]

B11[p_, m1_, m2_] := 1/18 - 2/3 B1[p, m1, m2] -
  2/3 (m1 - m2) DeltaB1[p, m1, m2] - 1/3 m1 DeltaB0[p, m1, m2]

(*
B11[p_, m1_, m2_] := 1/(3 p) (
  (p - 3 (m1 + m2))/6 + A0[m2] - m1 B0[p, m1, m2] -
  2 (p + m1 - m2) B1[p, m1, m2] )
*)


Derivative[1, 0, 0][B0] = DB0;
Derivative[1, 0, 0][B1] = DB1;
Derivative[1, 0, 0][B00] = DB00;
Derivative[1, 0, 0][B11] = DB11

DB0[0, m_, m_] := 1/(6 m)

DB0[0, m1_, m2_] :=
  1/(m1 - m2)^2 ((m1 + m2)/2 - A0[m2] + m2 B0[0, m1, m2]) /;
  ok[m1, m2]


DB1[0, m_, m_] = -1/(12 m)

DB1[0, m1_, m2_] =
  (2 m2 (B1[0, m1, m2] - B1[0, m2, m2])/(m1 - m2) - 1/3)/(m1 - m2) /;
  ok[m1, m2]

DB1[p_, m1_, m2_] := -DB0[p, m1, m2]/2 -
  (m2 - m1)/2 (DeltaB0[p, m1, m2] - DB0[p, m1, m2])/p



DB00[p_, m1_, m2_] :=
  1/6 (2 m1 DB0[p, m1, m2] + B1[p, m1, m2] +
    (p + m1 - m2) DB1[p, m1, m2] - 1/3)


DB11[p_, m1_, m2_] = D[B11[p, m1, m2], p] //Simplify


B0[p_, m1_, m2_] := B0[p, m2, m1] /; !OrderedQ[{m1, m2}];
DB0[p_, m1_, m2_] := DB0[p, m2, m1] /; !OrderedQ[{m1, m2}]

